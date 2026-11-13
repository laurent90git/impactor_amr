// Slip treatment at solid/fluid interfaces.
//
// Problem: in a diffuse-interface model the residual solid fraction gives the "fluid" cells a shear stiffness, so a
// solid sliding along a fluid drags it (no slip) and shear waves are sent into the fluid.
//
// A slip policy decides, cell by cell, whether the solid phases carry shear stress (chi = 1) or not (chi = 0), and
// what to do with their cobasis where they do not. Interface of a policy:
//   double chi(double alpha_s) const;                  // per solid phase s: 1 solid behaviour, 0 fluid behaviour
//   void reset(Mat3& A) const;                         // applied to each solid cobasis where chi = 0
//
// With chi = 0 the shear stress, the shear energy and the shear impedance of the cell vanish. Combined with the shear
// Riemann solver (riemann.hpp), whose traction is weighted by the shear impedances, a face between a solid cell and a
// fluid cell gets zero tangential traction: free slip, with a single (conservative) flux per face.
#pragma once

#include <cmath>

#include "../core/linalg.hpp"

namespace hvi
{
    /// Solid strength switched off below a threshold of total solid volume fraction; the cobasis is reset to a pure
    /// dilatation (the elastic energy it held becomes heat through the energy reset of the relaxation step).
    struct Fluidisation
    {
        double alpha_lim = 1e-2;

        [[nodiscard]] double chi(double alpha_solid) const // volume fraction of the solid phase
        {
            return alpha_solid >= alpha_lim ? 1. : 0.;
        }

        void reset(Mat3& A) const
        {
            const double J = std::cbrt(std::abs(det(A)));
            A              = J * identity3();
        }
    };

    /// No slip treatment: the solid keeps its strength whatever its volume fraction (the papers' raw model).
    struct NoSlipTreatment
    {
        [[nodiscard]] double chi(double) const
        {
            return 1.;
        }

        void reset(Mat3&) const
        {
        }
    };
} // namespace hvi
