// Driver: mesh, fields, initial condition, time loop (explicit SSP-RK2 + cell relaxation), multiresolution
// adaptation, outputs. Everything specific to the physics is in src/physics, to the numerics in src/numerics; this
// file only orchestrates them on a samurai mesh.
//
// Time step (rk = 2):   U1 = R(U - dt L(U)),   U <- R((U + U1 - dt L(U1)) / 2)
// with L the hyperbolic operator (flux_operator.hpp) and R the cell relaxation (physics/relaxation.hpp).
#pragma once

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

#include "../cases/cases.hpp"
#include "../io/npy.hpp"
#include "../physics/relaxation.hpp"
#include "boundary.hpp"
#include "flux_operator.hpp"

namespace hvi
{
    struct RunOptions
    {
        std::string out_dir = "out";
        std::string tag;               ///< suffix of the output files
        std::size_t level     = 0;     ///< finest level (0: the case's default)
        std::size_t min_level = 0;     ///< coarsest level (0 or = level: uniform mesh)
        double cfl            = 0.4;
        double tf             = -1.;   ///< final time (< 0: the case's)
        std::size_t nout      = 1;     ///< number of evenly spaced outputs (ignored if --paper-times)
        bool paper_times      = false; ///< output at the case's reference times
        std::string limiter   = "minmod";
        bool muscl_A          = true;
        int rk                = 2;
        double alim           = -1.;   ///< fluidisation threshold (< 0: the case's; 0: no slip treatment)
        int plastic           = -1;    ///< -1: the case's, 0: off, 1: on
        double eps            = -1.;   ///< residual volume fraction (< 0: the case's)
        double mr_eps         = 1e-3;
        double mr_reg         = 1.;
        std::size_t axis      = 0;     ///< direction of the 1D cases
        std::size_t log_every = 100;
        std::size_t max_steps = 100000000;
    };

    // ---------------------------------------------------------------------------------------------------------------
    // cell-wise helpers
    // ---------------------------------------------------------------------------------------------------------------
    template <class Field, class Model>
    typename Model::cons_t get_cons(const Field& U, const auto& cell)
    {
        typename Model::cons_t c;
        for (std::size_t q = 0; q < Model::nc; ++q)
        {
            c[q] = U[cell][q];
        }
        return c;
    }

    /// initial cell averages by sub-sampling of the case's point description
    template <std::size_t dim, class Model, class Field>
    void initialise(Field& U, const Model& model, const CaseDef<Model::N>& cd, double eps)
    {
        constexpr std::size_t N = Model::N, NS = Model::NS;
        const int ns            = dim == 1 ? 16 : (dim == 2 ? 6 : 4);
        samurai::for_each_cell(U.mesh(),
                               [&](auto& cell)
                               {
                                   std::array<double, N> frac{};
                                   Vec3 mom{};
                                   double mass = 0., p = 0., wsum = 0.;
                                   std::array<int, 3> k{0, 0, 0};
                                   const int ntot = static_cast<int>(std::pow(ns, dim));
                                   for (int s = 0; s < ntot; ++s)
                                   {
                                       int r = s;
                                       Vec3 x{};
                                       for (std::size_t d = 0; d < dim; ++d)
                                       {
                                           k[d] = r % ns;
                                           r /= ns;
                                           x[d] = cell.center(d) + cell.length * ((k[d] + 0.5) / ns - 0.5);
                                       }
                                       const PointState ps = cd.point(x);
                                       frac[ps.phase] += 1.;
                                       mass += cd.rho_init[ps.phase];
                                       for (int i = 0; i < 3; ++i)
                                       {
                                           mom[i] += cd.rho_init[ps.phase] * ps.v[i];
                                       }
                                       p += ps.p;
                                       wsum += 1.;
                                   }
                                   typename Model::prim_t w;
                                   for (std::size_t l = 0; l < N; ++l)
                                   {
                                       w.alpha[l] = eps + (1. - N * eps) * frac[l] / wsum;
                                       w.rho[l]   = cd.rho_init[l];
                                       w.p[l]     = p / wsum;
                                   }
                                   for (int i = 0; i < 3; ++i)
                                   {
                                       w.v[i] = mom[i] / mass;
                                   }
                                   for (std::size_t s = 0; s < NS; ++s)
                                   {
                                       // unstrained: g = I whatever the density; A = (rho/rho0)^(1/3) I
                                       w.A[s] = std::cbrt(w.rho[s] / model.mat[s].rho0) * identity3();
                                   }
                                   w.chi        = model.chi_of(w.alpha);
                                   const auto c = model.to_cons(w);
                                   for (std::size_t q = 0; q < Model::nc; ++q)
                                   {
                                       U[cell][q] = c[q];
                                   }
                               });
    }

    template <class Model, class Relaxation, class Field>
    RelaxationStats relax_all(Field& U, const Model& model, const Relaxation& relax)
    {
        RelaxationStats st;
        samurai::for_each_cell(U.mesh(),
                               [&](auto& cell)
                               {
                                   auto c = get_cons<Field, Model>(U, cell);
                                   relax(model, c, st);
                                   for (std::size_t q = 0; q < Model::nc; ++q)
                                   {
                                       U[cell][q] = c[q];
                                   }
                               });
        U.ghosts_updated() = false;
        return st;
    }

    // ---------------------------------------------------------------------------------------------------------------
    // the run
    // ---------------------------------------------------------------------------------------------------------------
    template <std::size_t dim, class Model>
    int run_case(Model model, CaseDef<Model::N> cd, RunOptions opt)
    {
        constexpr std::size_t N = Model::N, NS = Model::NS, nc = Model::nc;
        using L                 = typename Model::layout_t;

        // ---- options resolved against the case defaults --------------------------------------------------------------
        const std::size_t level = opt.level > 0 ? opt.level : cd.level;
        const std::size_t lmin  = (opt.min_level > 0 && opt.min_level < level) ? opt.min_level : level;
        const double tf         = opt.tf > 0. ? opt.tf : cd.tf;
        const double eps        = opt.eps > 0. ? opt.eps : cd.eps;
        model.slip.alpha_lim    = opt.alim >= 0. ? opt.alim : cd.alim;
        CellRelaxation<Model> relaxation;
        relaxation.plastic = opt.plastic < 0 ? cd.plastic : opt.plastic == 1;
        relaxation.rho_floor = 1e300;
        for (std::size_t l = 0; l < N; ++l)
        {
            relaxation.rho_floor = std::min(relaxation.rho_floor, 1e-3 * cd.rho_init[l]);
        }
        MusclPrimitive<Model> recon;
        recon.limiter       = limiter_from_string(opt.limiter);
        recon.reconstruct_A = opt.muscl_A;

        std::vector<double> out_times;
        if (opt.paper_times && !cd.times.empty())
        {
            for (double t : cd.times)
            {
                if (t <= tf * (1. + 1e-12))
                {
                    out_times.push_back(t);
                }
            }
        }
        else
        {
            for (std::size_t k = 1; k <= opt.nout; ++k)
            {
                out_times.push_back(tf * static_cast<double>(k) / static_cast<double>(opt.nout));
            }
        }

        // ---- mesh -----------------------------------------------------------------------------------------------------
        Vec3 lo = cd.lo, hi = cd.hi;
        if (cd.one_d) // extent of the case (given along x) on the axis, one level-0 cell across the other directions
        {
            lo = {0., 0., 0.};
            hi = {cd.base, cd.base, cd.base};
            lo[opt.axis] = cd.lo[0];
            hi[opt.axis] = cd.hi[0];
        }
        xt::xtensor_fixed<double, xt::xshape<dim>> blo, bhi;
        for (std::size_t d = 0; d < dim; ++d)
        {
            blo[d] = lo[d];
            bhi[d] = hi[d];
        }
        auto config = samurai::mesh_config<dim>().min_level(lmin).max_level(level).max_stencil_radius(2).disable_args_parse();
        config.scaling_factor(cd.base);
        auto mesh  = samurai::mra::make_mesh(samurai::Box<double, dim>(blo, bhi), config);
        auto U     = samurai::make_vector_field<double, nc>("U", mesh);
        auto U1    = samurai::make_vector_field<double, nc>("U1", mesh);
        auto Unp1  = samurai::make_vector_field<double, nc>("Unp1", mesh);
        auto ind   = samurai::make_vector_field<double, 3 + dim>("indicator", mesh);
        samurai::make_bc<Transmissive>(U);
        U1.copy_bc_from(U);
        Unp1.copy_bc_from(U);
        using mesh_id_t = typename decltype(mesh)::mesh_id_t;

        // ---- operator ------------------------------------------------------------------------------------------------
        auto op = make_flux_operator<decltype(U)>(model, recon);

        // ---- adaptation (multiresolution on a normalised indicator) ---------------------------------------------------
        const bool adaptive = lmin < level;
        double rho_ref = 1., v_ref = 1., p_ref = 1.;
        for (std::size_t l = 0; l < N; ++l)
        {
            rho_ref = std::max(rho_ref, cd.rho_init[l]);
        }
        {
            double vmax = 0.;
            const int nprobe = 50;
            for (int i = 0; i <= nprobe; ++i) // rough probe of the initial velocities
            {
                for (int j = 0; j <= (dim > 1 ? nprobe : 0); ++j)
                {
                    Vec3 x{lo[0] + (hi[0] - lo[0]) * i / nprobe, lo[1] + (hi[1] - lo[1]) * j / nprobe, 0.5 * (lo[2] + hi[2])};
                    const auto ps = cd.point(x);
                    vmax          = std::max({vmax, std::abs(ps.v[0]), std::abs(ps.v[1]), std::abs(ps.v[2])});
                }
            }
            v_ref = std::max(vmax, 100.);
            p_ref = rho_ref * v_ref * 3000.;
        }
        double rho_min_init = rho_ref;
        for (std::size_t l = 0; l < N; ++l)
        {
            rho_min_init = std::min(rho_min_init, cd.rho_init[l]);
        }
        const double log_rho_range = std::max(std::log(rho_ref / rho_min_init), 1.);
        auto fill_indicator        = [&]()
        {
            samurai::for_each_cell(mesh,
                                   [&](auto& cell)
                                   {
                                       const auto w = model.to_prim(get_cons<decltype(U), Model>(U, cell));
                                       double rho   = 0.;
                                       for (std::size_t l = 0; l < N; ++l)
                                       {
                                           rho += w.alpha[l] * w.rho[l];
                                       }
                                       double ap = 0.;
                                       for (std::size_t l = 0; l < N; ++l)
                                       {
                                           ap += w.alpha[l] * w.p[l];
                                       }
                                       ind[cell][0] = model.solid_fraction(w.alpha);
                                       // log: the shocks in the air are seen next to the solid/air density jump
                                       ind[cell][1] = std::log(rho) / log_rho_range;
                                       ind[cell][2] = ap / p_ref;
                                       for (std::size_t d = 0; d < dim; ++d)
                                       {
                                           ind[cell][3 + d] = w.v[d] / v_ref;
                                       }
                                   });
        };
        auto mra_cfg  = samurai::mra_config().epsilon(opt.mr_eps).regularity(opt.mr_reg);
        auto mr_adapt = samurai::make_MRAdapt(ind);
        auto adapt    = [&]()
        {
            fill_indicator();
            mr_adapt(mra_cfg, U);
            U1.resize();
            Unp1.resize();
        };

        // ---- initial condition ----------------------------------------------------------------------------------------
        initialise<dim>(U, model, cd, eps);
        if (adaptive)
        {
            for (std::size_t k = lmin; k <= level; ++k)
            {
                adapt();
                initialise<dim>(U, model, cd, eps); // exact values on the adapted mesh
            }
        }
        relax_all(U, model, relaxation);

        // ---- outputs ----------------------------------------------------------------------------------------------------
        std::filesystem::create_directories(opt.out_dir);
        const std::string prefix = fmt::format("{}/{}{}", opt.out_dir, cd.name, opt.tag.empty() ? "" : "_" + opt.tag);
        std::ofstream log(prefix + ".log");
        log << "# step t dt cells mass_0..mass_N-1 momentum_0..2 energy floored yielded\n";
        {
            std::ofstream hdr(prefix + "_columns.txt");
            hdr << "x y z h level rho alpha_solid u v w p s11 s22 s12 a21 chi vonmises s33 s13 s23\n";
            hdr << fmt::format("# case {} dim {} levels {}..{} base {} limiter {} rk {} cfl {} alim {} plastic {} eps {} N {} NS {}\n", cd.name,
                               dim, lmin, level, cd.base, opt.limiter, opt.rk, opt.cfl, model.slip.alpha_lim, relaxation.plastic, eps, N, NS);
            for (std::size_t l = 0; l < N; ++l)
            {
                hdr << fmt::format("# phase {}: {} gamma {} pinf {} mu {} rho0 {} yield {}\n", l, model.mat[l].name, model.mat[l].gamma,
                                   model.mat[l].pinf, model.mat[l].mu, model.mat[l].rho0, model.mat[l].yield);
            }
        }
        std::size_t iout = 0;
        auto output      = [&](double t)
        {
            std::vector<double> a;
            constexpr std::size_t ncol = 20;
            a.reserve(mesh.nb_cells(mesh_id_t::cells) * ncol);
            samurai::for_each_cell(mesh,
                                   [&](auto& cell)
                                   {
                                       const auto w  = model.to_prim(get_cons<decltype(U), Model>(U, cell));
                                       const auto mx = model.mixture(w, 0);
                                       double ap     = 0.;
                                       for (std::size_t l = 0; l < N; ++l)
                                       {
                                           ap += w.alpha[l] * w.p[l];
                                       }
                                       Mat3 S = mx.sigma;
                                       for (int i = 0; i < 3; ++i)
                                       {
                                           S[i][i] += ap;
                                       }
                                       const double vm = std::sqrt(1.5 * ddot(S, S));
                                       const double a21 = NS > 0 ? w.A[0][1][0] : 0.;
                                       Vec3 x{};
                                       for (std::size_t d = 0; d < dim; ++d)
                                       {
                                           x[d] = cell.center(d);
                                       }
                                       for (double v : {x[0], x[1], x[2], static_cast<double>(cell.length), static_cast<double>(cell.level), mx.rho,
                                                        model.solid_fraction(w.alpha), w.v[0], w.v[1], w.v[2], ap, mx.sigma[0][0],
                                                        mx.sigma[1][1], mx.sigma[0][1], a21, NS > 0 ? w.chi[0] : 0., vm, mx.sigma[2][2], mx.sigma[0][2],
                                                        mx.sigma[1][2]})
                                       {
                                           a.push_back(v);
                                       }
                                   });
            write_npy(fmt::format("{}_{:03d}.npy", prefix, iout), a, {a.size() / ncol, ncol});
            std::cout << fmt::format("  output {} at t = {:.4e} s ({} cells)", iout, t, mesh.nb_cells(mesh_id_t::cells)) << std::endl;
            ++iout;
        };

        auto diagnostics = [&](std::size_t step, double t, double dt, const RelaxationStats& st)
        {
            std::array<double, N> mass{};
            Vec3 mom{};
            double energy = 0.;
            samurai::for_each_cell(mesh,
                                   [&](auto& cell)
                                   {
                                       const double vol = std::pow(cell.length, dim);
                                       for (std::size_t l = 0; l < N; ++l)
                                       {
                                           mass[l] += U[cell][L::arho(l)] * vol;
                                       }
                                       for (std::size_t i = 0; i < 3; ++i)
                                       {
                                           mom[i] += U[cell][L::mom(i)] * vol;
                                       }
                                       energy += U[cell][L::energy] * vol;
                                   });
            std::string s = fmt::format("{} {:.10e} {:.4e} {}", step, t, dt, mesh.nb_cells(mesh_id_t::cells));
            for (double m : mass)
            {
                s += fmt::format(" {:.12e}", m);
            }
            s += fmt::format(" {:.12e} {:.12e} {:.12e} {:.12e} {} {}", mom[0], mom[1], mom[2], energy, st.floored, st.yielded);
            log << s << "\n";
            return s;
        };

        // ---- time loop ------------------------------------------------------------------------------------------------
        const double h  = mesh.min_cell_length();
        double t        = 0.;
        std::size_t nt = 0, next_out = 0;
        const auto wall0 = std::chrono::steady_clock::now();
        std::cout << fmt::format("case {} ({}D, N = {}, NS = {}): levels {}..{} (h = {:.4g} m), {} cells, tf = {:.4e} s", cd.name, dim, N,
                                 NS, lmin, level, h, mesh.nb_cells(mesh_id_t::cells), tf)
                  << std::endl;
        if (!opt.paper_times && opt.nout > 1)
        {
            output(0.);
        }
        RelaxationStats stats;
        double smax0 = 1.;
        while (t < tf * (1. - 1e-12) && nt < opt.max_steps)
        {
            if (adaptive)
            {
                adapt();
                const auto s = relax_all(U, model, relaxation);
                stats.floored += s.floored;
                stats.voided += s.voided;
                stats.vacuum += s.vacuum;
                stats.failed += s.failed;
            }
            double smax = 0.;
            bool finite = true;
            typename Model::prim_t wmax;
            Vec3 xmax{};
            samurai::for_each_cell(mesh,
                                   [&](auto& cell)
                                   {
                                       const auto w = model.to_prim(get_cons<decltype(U), Model>(U, cell));
                                       const double s = model.template cfl_speed<dim>(w);
                                       finite         = finite && std::isfinite(s);
                                       if (s > smax)
                                       {
                                           smax = s;
                                           wmax = w;
                                           for (std::size_t d = 0; d < dim; ++d)
                                           {
                                               xmax[d] = cell.center(d);
                                           }
                                       }
                                   });
            if (nt == 0)
            {
                smax0 = smax;
            }
            if (!finite || smax > 1e3 * smax0)
            {
                // the time step collapsed: report the cell that limits it
                std::cerr << fmt::format("{} at step {}, t = {:.6e}: max speed {:.4e} (initial {:.4e}) at ({:.5f}, {:.5f}, {:.5f})",
                                         finite ? "time step collapse" : "non-finite state", nt, t, smax, smax0, xmax[0], xmax[1], xmax[2])
                          << std::endl;
                for (std::size_t l = 0; l < N; ++l)
                {
                    std::cerr << fmt::format("  phase {}: alpha {:.6e} rho {:.6e} p {:.6e}", l, wmax.alpha[l], wmax.rho[l], wmax.p[l]) << std::endl;
                }
                std::cerr << fmt::format("  v ({:.4e}, {:.4e}, {:.4e}) chi {}", wmax.v[0], wmax.v[1], wmax.v[2], NS > 0 ? wmax.chi[0] : 0.) << std::endl;
                for (std::size_t s = 0; s < NS; ++s)
                {
                    const auto& A = wmax.A[s];
                    std::cerr << fmt::format("  A[{}] = [[{:.4e} {:.4e} {:.4e}] [{:.4e} {:.4e} {:.4e}] [{:.4e} {:.4e} {:.4e}]] det {:.4e}", s, A[0][0],
                                             A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2], det(A))
                              << std::endl;
                }
                output(t);
                return 1;
            }
            double dt = opt.cfl * h / smax;
            if (next_out < out_times.size())
            {
                dt = std::min(dt, out_times[next_out] - t);
            }
            dt = std::min(dt, tf - t);

            RelaxationStats st;
            if (opt.rk == 1)
            {
                Unp1 = U - dt * op(U);
                st   = relax_all(Unp1, model, relaxation);
            }
            else
            {
                U1      = U - dt * op(U);
                auto s1 = relax_all(U1, model, relaxation);
                Unp1    = 0.5 * (U + (U1 - dt * op(U1)));
                st      = relax_all(Unp1, model, relaxation);
                st.floored += s1.floored;
                st.yielded += s1.yielded;
                st.voided += s1.voided;
                st.vacuum += s1.vacuum;
                st.failed += s1.failed;
            }
            samurai::swap(U, Unp1);
            t += dt;
            ++nt;
            stats.floored += st.floored;
            stats.voided += st.voided;
            stats.vacuum += st.vacuum;
            stats.failed += st.failed;
            stats.yielded = st.yielded;
            const bool is_out = next_out < out_times.size() && t >= out_times[next_out] * (1. - 1e-12);
            if (nt % opt.log_every == 0 || is_out || nt == 1)
            {
                const auto s     = diagnostics(nt, t, dt, stats);
                const double wct = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
                std::cout << fmt::format("step {:7d} t = {:.4e} dt = {:.3e} cells {:8d} floored {:6d} voided {:7d} vacuum {:5d} failed {:6d} yielded {:6d}  ({:.0f} s)",
                                         nt, t, dt, mesh.nb_cells(mesh_id_t::cells), stats.floored, stats.voided, stats.vacuum, stats.failed,
                                         st.yielded, wct)
                          << std::endl;
                log.flush();
            }
            if (is_out)
            {
                output(t);
                ++next_out;
            }
        }
        const double wct = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
        std::cout << fmt::format("done: {} steps, t = {:.6e} s, {:.1f} s wall", nt, t, wct) << std::endl;
        return 0;
    }
} // namespace hvi
