// Multi-material hyperelastic diffuse-interface model of Ndanou, Favrie & Gavrilyuk (JCP 295, 2015), in its
// pressure non-equilibrium form (one velocity, N phases, instantaneous pressure relaxation).
//
// Phases 0 .. NS-1 are solids (they carry a cobasis A = F^{-1}), phases NS .. N-1 are fluids. The state is always
// three-dimensional (3 velocity components, 3x3 cobasis) whatever the dimension of the mesh.
//
// Conservative vector (Layout<N, NS>):
//   alpha_l                 N   volume fractions                     (non-conservative)
//   alpha_l rho_l           N   partial densities                     (conservative)
//   rho v                   3   momentum                              (conservative)
//   rho E                   1   total energy                          (conservative)
//   alpha_l rho_l e^h_l     N   hydrodynamic energies                 (non-conservative)
//   A_s                     9NS cobasis of each solid, row-major,  A[b][i] = d X^b / d x_i   (non-conservative)
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "../core/linalg.hpp"
#include "material.hpp"
#include "shear_energy.hpp"
#include "slip.hpp"

namespace hvi
{
    template <std::size_t N_, std::size_t NS_>
    struct Layout
    {
        static constexpr std::size_t N  = N_;
        static constexpr std::size_t NS = NS_;
        static_assert(NS <= N && N >= 1);

        static constexpr std::size_t alpha(std::size_t l) { return l; }
        static constexpr std::size_t arho(std::size_t l) { return N + l; }
        static constexpr std::size_t mom(std::size_t i) { return 2 * N + i; }
        static constexpr std::size_t energy = 2 * N + 3;
        static constexpr std::size_t arhoe(std::size_t l) { return 2 * N + 4 + l; }
        static constexpr std::size_t A(std::size_t s, std::size_t b, std::size_t i) { return 3 * N + 4 + 9 * s + 3 * b + i; }
        static constexpr std::size_t size = 3 * N + 4 + 9 * NS;
    };

    /// Primitive state of a cell or of a reconstructed face state.
    template <std::size_t N, std::size_t NS>
    struct Prim
    {
        std::array<double, N> alpha{}, rho{}, p{};
        Vec3 v{};
        std::array<Mat3, NS> A{};
        std::array<double, NS> chi{}; ///< strength switch of each solid (slip policy), never reconstructed
    };

    /// Mixture quantities derived from a primitive state, for one direction d.
    template <std::size_t N, std::size_t NS>
    struct Mixture
    {
        double rho = 0.;  ///< mixture density
        Mat3 sigma{};     ///< Cauchy stress
        double rhoE = 0.; ///< total energy per unit volume
        double cl = 0.;   ///< longitudinal wave speed (upper bound) in direction d
        double ct = 0.;   ///< transverse wave speed in direction d (0 where chi = 0)
    };

    template <std::size_t N_, std::size_t NS_, class ShearEnergy_ = NeoHookeanSeparable, class Slip_ = Fluidisation>
    class Model
    {
      public:

        static constexpr std::size_t N  = N_;
        static constexpr std::size_t NS = NS_;
        using layout_t                  = Layout<N, NS>;
        using ShearEnergy               = ShearEnergy_;
        using Slip                      = Slip_;
        using prim_t                    = Prim<N, NS>;
        using mixture_t                 = Mixture<N, NS>;
        using cons_t                    = std::array<double, layout_t::size>;
        static constexpr std::size_t nc = layout_t::size;

        std::array<Material, N> mat;
        Slip slip;
        double alpha_min = 1e-12; ///< floor on the volume fractions (robustness of predicted / reconstructed states)

        explicit Model(const std::array<Material, N>& materials, Slip s = {})
            : mat(materials)
            , slip(s)
        {
            for (std::size_t l = 0; l < N; ++l)
            {
                // a phase with a cobasis may have mu = 0 (fluid limit of a solid), a phase without one may not have mu > 0
                if (l >= NS && mat[l].is_solid())
                {
                    throw std::invalid_argument("Model: the NS first phases must be the solids, the others fluids (" + mat[l].name + ")");
                }
            }
        }

        // -------------------------------------------------------------------------------------------------------------
        // conversions
        // -------------------------------------------------------------------------------------------------------------
        /// strength switch of each solid phase: a solid carries shear stress only where its own volume fraction is
        /// above the threshold of the slip policy (residual solids inside other materials have no strength)
        [[nodiscard]] std::array<double, NS> chi_of(const std::array<double, N>& alpha) const
        {
            std::array<double, NS> c{};
            for (std::size_t s = 0; s < NS; ++s)
            {
                c[s] = slip.chi(alpha[s]);
            }
            return c;
        }

        [[nodiscard]] double solid_fraction(const std::array<double, N>& alpha) const
        {
            double a = 0.;
            for (std::size_t s = 0; s < NS; ++s)
            {
                a += alpha[s];
            }
            return a;
        }

        /// conservative -> primitive; robust to slightly inconsistent states (predicted ghosts): the volume fractions
        /// are clipped and normalised, the phase pressures kept above -p_inf.
        template <class C>
        [[nodiscard]] prim_t to_prim(const C& u) const
        {
            prim_t w;
            double sa = 0.;
            for (std::size_t l = 0; l < N; ++l)
            {
                w.alpha[l] = std::clamp(static_cast<double>(u[layout_t::alpha(l)]), alpha_min, 1.);
                sa += w.alpha[l];
            }
            double rho = 0.;
            for (std::size_t l = 0; l < N; ++l)
            {
                w.alpha[l] /= sa;
                const double ar = std::max(static_cast<double>(u[layout_t::arho(l)]), 1e-30);
                rho += ar;
                w.rho[l] = ar / w.alpha[l];
                w.p[l]   = mat[l].pressure(static_cast<double>(u[layout_t::arhoe(l)]) / w.alpha[l]);
                w.p[l]   = std::max(w.p[l], mat[l].pmin() + 1e-9 * (1. + mat[l].pinf));
            }
            for (std::size_t i = 0; i < 3; ++i)
            {
                w.v[i] = u[layout_t::mom(i)] / rho;
            }
            for (std::size_t s = 0; s < NS; ++s)
            {
                for (std::size_t b = 0; b < 3; ++b)
                {
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        w.A[s][b][i] = u[layout_t::A(s, b, i)];
                    }
                }
            }
            w.chi = chi_of(w.alpha);
            return w;
        }

        [[nodiscard]] cons_t to_cons(const prim_t& w) const
        {
            cons_t u{};
            const auto mx = mixture(w, 0);
            double rho    = 0.;
            for (std::size_t l = 0; l < N; ++l)
            {
                u[layout_t::alpha(l)] = w.alpha[l];
                u[layout_t::arho(l)]  = w.alpha[l] * w.rho[l];
                u[layout_t::arhoe(l)] = w.alpha[l] * mat[l].rho_eh(w.p[l]);
                rho += w.alpha[l] * w.rho[l];
            }
            for (std::size_t i = 0; i < 3; ++i)
            {
                u[layout_t::mom(i)] = rho * w.v[i];
            }
            u[layout_t::energy] = mx.rhoE;
            for (std::size_t s = 0; s < NS; ++s)
            {
                for (std::size_t b = 0; b < 3; ++b)
                {
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        u[layout_t::A(s, b, i)] = w.A[s][b][i];
                    }
                }
            }
            return u;
        }

        // -------------------------------------------------------------------------------------------------------------
        // mixture quantities
        // -------------------------------------------------------------------------------------------------------------
        /// specific shear energy of solid s (0 where chi = 0)
        [[nodiscard]] double shear_energy(const prim_t& w, std::size_t s) const
        {
            if (w.chi[s] == 0.)
            {
                return 0.;
            }
            return ShearEnergy::energy(ShearEnergy::evaluate(w.A[s]), mat[s]);
        }

        [[nodiscard]] mixture_t mixture(const prim_t& w, std::size_t d) const
        {
            mixture_t m;
            double ap = 0., cf2 = 0., kin = 0., reh = 0., ree = 0.;
            for (std::size_t l = 0; l < N; ++l)
            {
                m.rho += w.alpha[l] * w.rho[l];
            }
            Mat3 K{}; // acoustic tensor per unit mass
            for (std::size_t l = 0; l < N; ++l)
            {
                const double Y = w.alpha[l] * w.rho[l] / m.rho;
                ap += w.alpha[l] * w.p[l];
                cf2 += Y * mat[l].c2(w.rho[l], w.p[l]);
                reh += w.alpha[l] * mat[l].rho_eh(w.p[l]);
            }
            for (std::size_t i = 0; i < 3; ++i)
            {
                m.sigma[i][i] = -ap;
                kin += w.v[i] * w.v[i];
            }
            K[d][d] = cf2;
            double chimax = 0.;
            for (std::size_t s = 0; s < NS; ++s)
            {
                chimax = std::max(chimax, w.chi[s]);
                if (w.chi[s] > 0.)
                {
                    const double ar = w.alpha[s] * w.rho[s];
                    const auto es   = ShearEnergy::evaluate(w.A[s]);
                    m.sigma += (w.chi[s] * w.alpha[s]) * ShearEnergy::deviator(es, mat[s], w.rho[s]);
                    ree += w.chi[s] * ar * ShearEnergy::energy(es, mat[s]);
                    K += (w.chi[s] * ar / m.rho) * ShearEnergy::acoustic(es, mat[s], d);
                }
            }
            m.rhoE = reh + ree + 0.5 * m.rho * kin;
            // longitudinal speed: Gershgorin bound of the largest eigenvalue of K (safe for the HLL wave speeds)
            double lmax = 0.;
            for (int i = 0; i < 3; ++i)
            {
                lmax = std::max(lmax, std::abs(K[i][0]) + std::abs(K[i][1]) + std::abs(K[i][2]));
            }
            m.cl = std::sqrt(lmax);
            if (chimax > 0.)
            {
                const std::size_t t1 = (d + 1) % 3, t2 = (d + 2) % 3;
                m.ct                 = std::sqrt(std::max(std::max(K[t1][t1], K[t2][t2]) + std::abs(K[t1][t2]), 0.));
            }
            return m;
        }

        /// largest |v_d| + c_l over the directions d < dim, summed (unsplit CFL condition)
        template <std::size_t dim>
        [[nodiscard]] double cfl_speed(const prim_t& w) const
        {
            double s = 0.;
            for (std::size_t d = 0; d < dim; ++d)
            {
                s += std::abs(w.v[d]) + mixture(w, d).cl;
            }
            return s;
        }
    };
} // namespace hvi
