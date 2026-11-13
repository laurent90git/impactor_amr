// Cell-local relaxation step, applied after each stage of the hyperbolic update (and after mesh adaptation):
//   1. sanitise: volume fractions clipped and normalised; vacuum guard (density floor, velocity cap; counted in
//      `vacuum`: these cells do not conserve mass and energy);
//   2. slip policy: where chi = 0, each solid cobasis is reset (slip.hpp);
//   3. plasticity: each solid cobasis relaxed to the yield surface (plasticity.hpp);
//   4. instantaneous pressure relaxation (Saurel et al. 2009 sec. 3.3, Ndanou 2015 sec. 5.3): solve for p
//        f(p) = sum_l alpha_l (p0_l + gamma_l pinf_l + (gamma_l - 1) p) / (gamma_l (p + pinf_l)) - 1 = 0
//      and update the volume fractions;
//   5. energy reset from the total energy (mixture EOS, Ndanou 2015 sec. 5.4). Where the energy-based pressure is
//      not admissible (strong expansion: the stiff phases were drained), a void is opened: at the floor pressure the
//      volume of the stiff phases is reduced and the gas fills the rest, so that the total energy is kept exactly
//      (counted in `voided`). Only if that is impossible (no gas phase) is the relaxed pressure kept and the total
//      energy changed (counted in `floored`).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "model.hpp"
#include "plasticity.hpp"

namespace hvi
{
    struct RelaxationStats
    {
        std::size_t floored = 0; ///< cells where the relaxed pressure was kept (energy not conserved locally)
        std::size_t yielded = 0; ///< solid cobases relaxed to the yield surface
        std::size_t voided  = 0; ///< cells where a void was opened to keep the energy (see step 5)
        std::size_t vacuum  = 0; ///< cells refilled by the vacuum guard (step 1b; mass and energy not conserved)
    };

    template <class Model, class Plasticity = VonMisesInstantaneous<typename Model::ShearEnergy>>
    struct CellRelaxation
    {
        using L      = typename Model::layout_t;
        using cons_t = typename Model::cons_t;

        Plasticity plasticity{};
        bool plastic      = true; ///< run-time switch of the plasticity step
        bool void_opening = true; ///< step 5: open a void rather than inject energy where the pressure is inadmissible
        double rho_floor  = 1e-3; ///< step 1b: density floor [kg/m3]
        double v_cap      = 2e4;  ///< step 1b: velocity cap [m/s]

        /// u: the conservative state of one cell (any indexable type)
        template <class C>
        void operator()(const Model& model, C& u, RelaxationStats& stats) const
        {
            constexpr std::size_t N = Model::N, NS = Model::NS;
            const auto& mat         = model.mat;

            // 1. sanitise
            std::array<double, N> al{};
            double sa = 0.;
            for (std::size_t l = 0; l < N; ++l)
            {
                al[l] = std::clamp(static_cast<double>(u[L::alpha(l)]), model.alpha_min, 1.);
                sa += al[l];
                u[L::arho(l)] = std::max(static_cast<double>(u[L::arho(l)]), 1e-30);
            }
            for (auto& a : al)
            {
                a /= sa;
            }

            // 1b. vacuum guard: the scheme is not positivity-preserving; a cell drained to (near) vacuum is refilled
            // with its lightest phase up to rho_floor, its velocity capped, its pressure set to the floor pressure.
            {
                double rho = 0.;
                for (std::size_t l = 0; l < N; ++l)
                {
                    rho += u[L::arho(l)];
                }
                double v2 = 0.;
                for (std::size_t i = 0; i < 3; ++i)
                {
                    v2 += sq(u[L::mom(i)] / rho);
                }
                if (rho < rho_floor || v2 > v_cap * v_cap)
                {
                    std::size_t lg = 0; // lightest phase: lowest p_inf
                    for (std::size_t l = 1; l < N; ++l)
                    {
                        if (mat[l].pinf < mat[lg].pinf)
                        {
                            lg = l;
                        }
                    }
                    if (rho < rho_floor)
                    {
                        u[L::arho(lg)] += rho_floor - rho;
                        rho = rho_floor;
                    }
                    const double vn = std::sqrt(v2);
                    const double f  = vn > v_cap ? v_cap / vn : 1.;
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        u[L::mom(i)] *= f;
                    }
                    // floor pressure in every phase; the total energy is recomputed in step 5 (bad branch)
                    double pfl = -1e300;
                    for (std::size_t l = 0; l < N; ++l)
                    {
                        pfl = std::max(pfl, mat[l].pmin());
                    }
                    double reh = 0., kin = 0.;
                    for (std::size_t l = 0; l < N; ++l)
                    {
                        u[L::arhoe(l)] = al[l] * mat[l].rho_eh(pfl + 1e-9 * (1. + std::abs(pfl)));
                        reh += u[L::arhoe(l)];
                    }
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        kin += 0.5 * sq(u[L::mom(i)]) / rho;
                    }
                    u[L::energy] = reh + kin; // shear energy added back in step 5 if the solid has strength
                    ++stats.vacuum;
                }
            }

            // 2 - 3. slip policy, plasticity
            auto get_A = [&](std::size_t s)
            {
                Mat3 A;
                for (std::size_t b = 0; b < 3; ++b)
                {
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        A[b][i] = u[L::A(s, b, i)];
                    }
                }
                return A;
            };
            auto set_A = [&](std::size_t s, const Mat3& A)
            {
                for (std::size_t b = 0; b < 3; ++b)
                {
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        u[L::A(s, b, i)] = A[b][i];
                    }
                }
            };
            const auto chi0 = model.chi_of(al);
            for (std::size_t s = 0; s < NS; ++s)
            {
                Mat3 A = get_A(s);
                if (chi0[s] == 0.)
                {
                    model.slip.reset(A);
                }
                else if (plastic && plasticity(A, mat[s], u[L::arho(s)] / al[s]))
                {
                    ++stats.yielded;
                }
                set_A(s, A);
            }

            // 4. pressure relaxation
            std::array<double, N> p0{};
            double pmin = -1e300, p = 0.;
            for (std::size_t l = 0; l < N; ++l)
            {
                pmin  = std::max(pmin, mat[l].pmin());
                p0[l] = std::max(mat[l].pressure(u[L::arhoe(l)] / al[l]), mat[l].pmin() + 1e-9 * (1. + mat[l].pinf));
                p += al[l] * p0[l];
            }
            const double ptol = 1e-9 * (1. + std::abs(pmin));
            if (N > 1)
            {
                p = std::max(p, pmin + ptol);
                for (int it = 0; it < 200; ++it)
                {
                    double f = -1., df = 0.;
                    for (std::size_t l = 0; l < N; ++l)
                    {
                        const double g = mat[l].gamma, pi = mat[l].pinf;
                        f += al[l] * (p0[l] + g * pi + (g - 1.) * p) / (g * (p + pi));
                        df -= al[l] * (p0[l] + pi) / (g * (p + pi) * (p + pi));
                    }
                    double pn = p - f / df;
                    if (pn <= pmin)
                    {
                        pn = 0.5 * (p + pmin);
                    }
                    const bool done = std::abs(pn - p) <= 1e-12 * (std::abs(p) + 1.);
                    p               = pn;
                    if (done)
                    {
                        break;
                    }
                }
                double s2 = 0.;
                for (std::size_t l = 0; l < N; ++l)
                {
                    const double g = mat[l].gamma, pi = mat[l].pinf;
                    al[l] *= (p0[l] + g * pi + (g - 1.) * p) / (g * (p + pi));
                    al[l] = std::max(al[l], model.alpha_min);
                    s2 += al[l];
                }
                for (auto& a : al)
                {
                    a /= s2;
                }
            }
            for (std::size_t l = 0; l < N; ++l)
            {
                u[L::alpha(l)] = al[l];
            }

            // 5. energy reset
            double rho = 0., kin = 0., ree = 0., a_g = 0., a_gp = 0.;
            for (std::size_t l = 0; l < N; ++l)
            {
                rho += u[L::arho(l)];
                const double g = mat[l].gamma;
                a_g += al[l] / (g - 1.);
                a_gp += al[l] * g * mat[l].pinf / (g - 1.);
            }
            for (std::size_t i = 0; i < 3; ++i)
            {
                kin += 0.5 * u[L::mom(i)] * u[L::mom(i)] / rho;
            }
            const auto chi = model.chi_of(al);
            for (std::size_t s = 0; s < NS; ++s)
            {
                using SE = typename Model::ShearEnergy;
                if (chi[s] > 0.)
                {
                    ree += u[L::arho(s)] * SE::energy(SE::evaluate(get_A(s)), mat[s]);
                }
            }
            const double reh_target = u[L::energy] - kin - ree; // hydrodynamic energy available
            double pe               = (reh_target - a_gp) / a_g;
            bool bad                = !(pe > pmin + ptol) || !std::isfinite(pe);
            if (bad && void_opening)
            {
                // Not enough energy for an admissible pressure with these volume fractions: the cell is in strong
                // expansion and its stiff phases (p_inf > pmin) were drained. Open the void: at the floor pressure
                // p* the stiff phases are compacted (volume fractions scaled by k < 1) and the gaseous ones fill the
                // rest, so that the mixture EOS matches the total energy exactly.
                const double ps = pmin + ptol;
                double as = 0., cs = 0., ag = 0., cg = 0.;
                for (std::size_t l = 0; l < N; ++l)
                {
                    const double c = al[l] * mat[l].rho_eh(ps);
                    if (mat[l].pinf > pmin)
                    {
                        as += al[l];
                        cs += c;
                    }
                    else
                    {
                        ag += al[l];
                        cg += c;
                    }
                }
                if (as > 0. && ag > 0.)
                {
                    const double Cg = cg / ag; // energy per unit volume of the gaseous phases at p*
                    const double k  = (reh_target - Cg) / (cs - as * Cg);
                    if (std::isfinite(k) && k > 0. && k < 1.)
                    {
                        const double m = (1. - k * as) / ag;
                        for (std::size_t l = 0; l < N; ++l)
                        {
                            al[l] *= mat[l].pinf > pmin ? k : m;
                            al[l] = std::max(al[l], model.alpha_min);
                            u[L::alpha(l)] = al[l];
                        }
                        pe  = ps;
                        bad = false;
                        ++stats.voided;
                    }
                }
            }
            const double pf = bad ? p : pe;
            double reh      = 0.;
            for (std::size_t l = 0; l < N; ++l)
            {
                u[L::arhoe(l)] = al[l] * mat[l].rho_eh(pf);
                reh += u[L::arhoe(l)];
            }
            if (bad)
            {
                u[L::energy] = reh + ree + kin;
                ++stats.floored;
            }
        }
    };
} // namespace hvi
