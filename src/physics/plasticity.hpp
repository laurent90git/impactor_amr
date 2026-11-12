// Plastic relaxation of a solid cobasis to the von Mises yield surface (Maxwell-type visco-plasticity in the limit of
// an instantaneous relaxation time), after Favrie & Gavrilyuk (2012) sec. 3.4 / 5.3 and Ndanou et al. (2015) sec. 5.2.
//
// Interface of a plasticity model (to bring another one, write a struct with the same operator()):
//   bool operator()(Mat3& A, const Material& m, double rho_l) const;   // relax A in place, true if it yielded
//
// VonMisesInstantaneous: with G = A^T A = U diag(lambda) U^T and b_a = lambda_a / (lambda_1 lambda_2 lambda_3)^(1/3),
// the principal deviatoric stresses S_a(b) are those of the shear energy. If S:S > 2/3 Y^2, the log-strains
// x_a = ln b_a follow the relaxation path dx_a/ds = S_a / |S| (Godunov's path: U and the left singular vectors are
// frozen, sum x_a = 0 so the density is unchanged, S:S decreases along it) until S:S = 2/3 Y^2. Then
//   A <- A U diag(exp(dx_a / 2)) U^T.
// The elastic energy released becomes heat through the energy reset (pressure_relaxation.hpp): no explicit source.
#pragma once

#include <cmath>

#include "../core/linalg.hpp"
#include "material.hpp"

namespace hvi
{
    template <class ShearEnergy>
    struct VonMisesInstantaneous
    {
        bool operator()(Mat3& A, const Material& m, double rho_l) const
        {
            if (!m.is_plastic())
            {
                return false;
            }
            const double limit = 2. / 3. * m.yield * m.yield;
            Vec3 lambda;
            Mat3 U;
            eigen_sym(ata(A), lambda, U);
            const double J3 = std::cbrt(std::abs(lambda[0] * lambda[1] * lambda[2]));
            Vec3 x0;
            for (int a = 0; a < 3; ++a)
            {
                x0[a] = std::log(std::max(lambda[a], 1e-300) / J3);
            }
            auto S_of = [&](const Vec3& x)
            { return ShearEnergy::principal_deviator({std::exp(x[0]), std::exp(x[1]), std::exp(x[2])}, m, rho_l); };
            auto f = [&](const Vec3& x)
            {
                const Vec3 S = S_of(x);
                return S[0] * S[0] + S[1] * S[1] + S[2] * S[2] - limit;
            };
            if (f(x0) <= 0.)
            {
                return false;
            }
            // march along the path (RK2, arc-length steps) until the yield surface is crossed, then bisect
            auto dir = [&](const Vec3& x)
            {
                const Vec3 S   = S_of(x);
                const double n = std::sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
                return Vec3{S[0] / n, S[1] / n, S[2] / n};
            };
            auto step = [&](const Vec3& x, double h)
            {
                const Vec3 k1 = dir(x);
                const Vec3 xm{x[0] + 0.5 * h * k1[0], x[1] + 0.5 * h * k1[1], x[2] + 0.5 * h * k1[2]};
                const Vec3 k2 = dir(xm);
                return Vec3{x[0] + h * k2[0], x[1] + h * k2[1], x[2] + h * k2[2]};
            };
            const double eps_y = m.yield / m.mu; // order of magnitude of the elastic strain at yield
            Vec3 x             = x0;
            for (int it = 0; it < 100000; ++it)
            {
                const double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
                const double h = std::max(0.1 * (r - 0.5 * eps_y), 0.05 * eps_y);
                const Vec3 xn  = step(x, h);
                if (f(xn) > 0.)
                {
                    x = xn;
                    continue;
                }
                double lo = 0., hi = h; // f(step(x, lo)) > 0 >= f(step(x, hi))
                for (int k = 0; k < 60 && hi - lo > 1e-12 * h; ++k)
                {
                    const double mid = 0.5 * (lo + hi);
                    (f(step(x, mid)) > 0. ? lo : hi) = mid;
                }
                x = step(x, hi);
                break;
            }
            Mat3 D{};
            for (int a = 0; a < 3; ++a)
            {
                D[a][a] = std::exp(0.5 * (x[a] - x0[a]));
            }
            A = A * (U * D * transpose(U));
            return true;
        }
    };

    /// No plasticity (purely elastic solids).
    struct NoPlasticity
    {
        bool operator()(Mat3&, const Material&, double) const
        {
            return false;
        }
    };
} // namespace hvi
