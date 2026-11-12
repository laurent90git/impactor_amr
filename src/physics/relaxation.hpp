// Cell-local relaxation step, applied after each stage of the hyperbolic update (and after mesh adaptation):
//   1. sanitise: volume fractions clipped and normalised;
//   2. slip policy: where chi = 0, each solid cobasis is reset (slip.hpp);
//   3. plasticity: each solid cobasis relaxed to the yield surface (plasticity.hpp);
//   4. instantaneous pressure relaxation (Saurel et al. 2009 sec. 3.3, Ndanou 2015 sec. 5.3): solve for p
//        f(p) = sum_l alpha_l (p0_l + gamma_l pinf_l + (gamma_l - 1) p) / (gamma_l (p + pinf_l)) - 1 = 0
//      and update the volume fractions;
//   5. energy reset from the total energy (mixture EOS, Ndanou 2015 sec. 5.4). Where the energy-based pressure is
//      not admissible (gas close to vacuum in an opening void), the relaxed pressure is kept and the total energy is
//      made consistent with it (counted in `floored`).
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
    };

    template <class Model, class Plasticity = VonMisesInstantaneous<typename Model::ShearEnergy>>
    struct CellRelaxation
    {
        using L      = typename Model::layout_t;
        using cons_t = typename Model::cons_t;

        Plasticity plasticity{};
        bool plastic = true; ///< run-time switch of the plasticity step

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
            const double chi0 = model.slip.chi(model.solid_fraction(al));
            for (std::size_t s = 0; s < NS; ++s)
            {
                Mat3 A = get_A(s);
                if (chi0 == 0.)
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
            const double chi = model.slip.chi(model.solid_fraction(al));
            if (chi > 0.)
            {
                for (std::size_t s = 0; s < NS; ++s)
                {
                    using SE = typename Model::ShearEnergy;
                    ree += u[L::arho(s)] * SE::energy(SE::evaluate(get_A(s)), mat[s]);
                }
            }
            const double pe = (u[L::energy] - kin - ree - a_gp) / a_g;
            const bool bad  = !(pe > pmin + ptol) || !std::isfinite(pe);
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
