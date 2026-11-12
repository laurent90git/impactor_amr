// Approximate Riemann solver at a face of normal direction d.
//
// Interface of a face solver (to bring another one, write a struct with the same operator()):
//   FaceResult<Model> operator()(const Model&, const prim_t& wl, const prim_t& wr, std::size_t d) const;
// It returns
//   - flux:  the numerical flux of every component. For the non-conservative components (alpha_l, alpha_l rho_l e_l,
//            A) it is the "advective" part  u*_d phi*  only;
//   - vstar: the velocity at the face (normal component along d, tangential ones from the shear solver). The
//            non-conservative products of the cell update (flux_operator.hpp) use its differences across the cell.
//
// LongitudinalHllcShear: sum of
//   - an HLLC solver of the longitudinal sub-system (mass, normal momentum, energy; the tangential momentum is only
//     advected), with the wave speeds of the full mixture (Ndanou 2015 sec. 5.1.1, Favrie 2009 sec. 4.2), the phase
//     energies along their Hugoniot (Saurel 2009 eq. II.8);
//   - a two-wave linear solver of each shear sub-system with the shear impedances Z = rho c_t (after Favrie,
//     Gavrilyuk, Ndanou 2014, sub-systems 2-3): with Z = 0 on one side (fluid, or fluidised solid), the tangential
//     traction vanishes and the face slides.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

#include "../physics/model.hpp"

namespace hvi
{
    template <class Model>
    struct FaceResult
    {
        typename Model::cons_t flux{};
        Vec3 vstar{};
    };

    template <class Model>
    struct LongitudinalHllcShear
    {
        using prim_t = typename Model::prim_t;
        using L      = typename Model::layout_t;

        FaceResult<Model> operator()(const Model& model, const prim_t& wl, const prim_t& wr, std::size_t d) const
        {
            constexpr std::size_t N = Model::N, NS = Model::NS;
            FaceResult<Model> res;
            auto& F = res.flux;

            const auto ml = model.mixture(wl, d), mr = model.mixture(wr, d);
            const double uL = wl.v[d], uR = wr.v[d];
            const double rL = ml.rho, rR = mr.rho;
            const double sL = ml.sigma[d][d], sR = mr.sigma[d][d];
            const double SL = std::min(uL - ml.cl, uR - mr.cl);
            const double SR = std::max(uL + ml.cl, uR + mr.cl);
            const double SM = (sL - sR + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR));

            // ---- longitudinal HLLC: upwind side K of the contact -----------------------------------------------------
            const bool left       = SM >= 0.;
            const prim_t& w       = left ? wl : wr;
            const auto& m         = left ? ml : mr;
            const double S        = left ? SL : SR;
            const double u        = w.v[d];
            const bool supersonic = left ? SL >= 0. : SR <= 0.;
            const double fac      = supersonic ? 1. : (S - u) / (S - SM);
            const double ustar    = supersonic ? u : SM;
            const double rho = m.rho, sig = m.sigma[d][d], E = m.rhoE;

            for (std::size_t l = 0; l < N; ++l)
            {
                const double ar         = w.alpha[l] * w.rho[l];
                F[L::arho(l)]           = ar * fac * ustar; // = ar u + S (ar fac - ar)
                F[L::alpha(l)]          = ustar * w.alpha[l];
                const Material& mat     = model.mat[l];
                const double rs         = w.rho[l] * fac;
                const double pK         = w.p[l] + mat.pinf;
                const double den        = (mat.gamma - 1.) * rs - (mat.gamma + 1.) * w.rho[l];
                double ps               = 0.;
                const double hug        = den < 0. ? pK * ((mat.gamma - 1.) * w.rho[l] - (mat.gamma + 1.) * rs) / den : -1.;
                ps                      = (den < 0. && hug > 0.) ? hug : pK * std::pow(fac, mat.gamma); // Hugoniot, else isentrope
                ps -= mat.pinf;
                F[L::arhoe(l)] = ustar * w.alpha[l] * mat.rho_eh(ps);
            }
            for (std::size_t i = 0; i < 3; ++i)
            {
                const double Fi = rho * u * w.v[i] - (i == d ? sig : 0.);
                const double Ui = rho * w.v[i];
                const double Us = rho * fac * (i == d ? SM : w.v[i]);
                F[L::mom(i)]    = supersonic ? Fi : Fi + S * (Us - Ui);
            }
            {
                const double FE = (E - sig) * u;
                const double Es = supersonic ? E : fac * (E + (SM - u) * (rho * SM - sig / (S - u)));
                F[L::energy]    = supersonic ? FE : FE + S * (Es - E);
            }
            for (std::size_t s = 0; s < NS; ++s)
            {
                for (std::size_t b = 0; b < 3; ++b)
                {
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        // the normal column behaves like a density across the acoustic waves, the others are advected
                        F[L::A(s, b, i)] = ustar * w.A[s][b][i] * (i == d ? fac : 1.);
                    }
                }
            }
            res.vstar[d] = ustar;

            // ---- shear sub-systems: tangential traction and velocity at the face -------------------------------------
            const double ZL = rL * ml.ct, ZR = rR * mr.ct, Zs = ZL + ZR;
            for (std::size_t t = 0; t < 3; ++t)
            {
                if (t == d)
                {
                    continue;
                }
                const double tL = ml.sigma[t][d], tR = mr.sigma[t][d];
                const double vL = wl.v[t], vR = wr.v[t];
                double vs = 0.5 * (vL + vR), ts = 0.;
                if (Zs > 0.)
                {
                    vs = (ZL * vL + ZR * vR + tR - tL) / Zs;
                    ts = (ZR * tL + ZL * tR + ZL * ZR * (vR - vL)) / Zs;
                }
                F[L::mom(t)] -= ts;
                F[L::energy] -= ts * vs;
                res.vstar[t] = vs;
            }
            return res;
        }
    };
} // namespace hvi
