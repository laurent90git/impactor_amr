// Shear (isochoric) part of the energy of a solid phase, as a function of its cobasis A = F^{-1}
// (row b of A is the cobasis vector e^b = grad X^b).
//
// Interface of a shear-energy model (to plug in another one, write a struct with the same members):
//   struct ElasticState { Mat3 g, g2; double j2; ... };
//   static ElasticState evaluate(const Mat3& A);
//   static double energy(const ElasticState&, const Material&);           // specific shear energy e^e [J/kg]
//   static Mat3 deviator(const ElasticState&, const Material&, double rho_l); // deviatoric stress S_l [Pa]
//   static Mat3 acoustic(const ElasticState&, const Material&, std::size_t d); // elastic acoustic tensor per unit mass
//
// NeoHookeanSeparable: Ndanou, Favrie, Gavrilyuk (2015) eq. (5):
//   G = A^T A,  g = G / |G|^(1/3),  j2 = tr(g^2),  e^e = mu / (8 rho0) (j2 - 3)
//   S = -2 rho de/dG G = -(mu rho_l / (2 rho0)) (g^2 - j2/3 I)
// Acoustic tensor (derived for this code, checked against a uniaxial finite-difference calculation):
//   D = mu/(2 rho0) [ g_dd g + (g e_d)(g e_d)^T + g^2 - ((g^2 e_d) e_d^T + e_d (g^2 e_d)^T)/3 + j2/9 e_d e_d^T ]
#pragma once

#include <cmath>
#include <cstddef>

#include "../core/linalg.hpp"
#include "material.hpp"

namespace hvi
{
    struct NeoHookeanSeparable
    {
        struct ElasticState
        {
            Mat3 g{}, g2{};
            double j2 = 3.;
        };

        static ElasticState evaluate(const Mat3& A)
        {
            ElasticState s;
            const Mat3 G     = ata(A);
            const double J23 = std::pow(std::abs(det(A)), 2. / 3.); // |G|^(1/3) = |det A|^(2/3)
            s.g              = (1. / J23) * G;
            s.g2             = s.g * s.g;
            s.j2             = trace(s.g2);
            return s;
        }

        static double energy(const ElasticState& s, const Material& m)
        {
            return m.mu / (8. * m.rho0) * (s.j2 - 3.);
        }

        static Mat3 deviator(const ElasticState& s, const Material& m, double rho_l)
        {
            const double k = -m.mu * rho_l / (2. * m.rho0);
            Mat3 S         = s.g2;
            for (int i = 0; i < 3; ++i)
            {
                S[i][i] -= s.j2 / 3.;
            }
            return k * S;
        }

        static Mat3 acoustic(const ElasticState& s, const Material& m, std::size_t d)
        {
            Mat3 D{};
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    D[i][j] = s.g[d][d] * s.g[i][j] + s.g[i][d] * s.g[j][d] + s.g2[i][j];
                }
            }
            for (int i = 0; i < 3; ++i)
            {
                D[i][d] -= s.g2[i][d] / 3.;
                D[d][i] -= s.g2[i][d] / 3.;
            }
            D[d][d] += s.j2 / 9.;
            return (m.mu / (2. * m.rho0)) * D;
        }

        /// principal deviatoric stresses from the eigenvalues b_a of g (product 1), used by the plastic relaxation
        static Vec3 principal_deviator(const Vec3& b, const Material& m, double rho_l)
        {
            const double k  = -m.mu * rho_l / (2. * m.rho0);
            const double j2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
            return {k * (b[0] * b[0] - j2 / 3.), k * (b[1] * b[1] - j2 / 3.), k * (b[2] * b[2] - j2 / 3.)};
        }
    };
} // namespace hvi
