// Test cases. A case gives the materials, the domain, and a point-wise description of the initial state: which phase
// fills the point, its velocity and pressure. Cells are initialised by sub-sampling (cell averages of the volume
// fractions, of the mass and of the momentum); a residual fraction eps of every other phase is added everywhere.
//
// 1D cases are written along the coordinate x[axis]; they run in any dimension (the other directions are uniform).
//
// Adding a case: add a branch in make_case_2phase / make_case_3phase and, if needed, a new model instantiation in
// app/main.cpp.
#pragma once

#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "../core/linalg.hpp"
#include "../physics/material.hpp"

namespace hvi
{
    struct PointState
    {
        std::size_t phase = 0; ///< material filling the point
        Vec3 v{};              ///< velocity
        double p = 1e5;        ///< pressure
    };

    template <std::size_t N>
    struct CaseDef
    {
        std::string name;
        std::array<Material, N> mat;
        std::array<double, N> rho_init{}; ///< initial phase densities
        Vec3 lo{0., 0., 0.}, hi{1., 1., 1.};
        double base       = 0.05; ///< cell size at level 0 (must divide the domain lengths)
        double eps        = 1e-6; ///< residual volume fraction of the absent phases
        double tf         = 1e-4;
        std::size_t level = 6;    ///< default finest level
        bool plastic      = false;
        double alim       = 1e-2;
        bool one_d        = false; ///< 1D case: only x[axis] matters
        std::vector<double> times; ///< reference output times of the paper (if any)
        std::function<PointState(const Vec3&)> point;
    };

    // -----------------------------------------------------------------------------------------------------------------
    // cases with one solid (phase 0) and air (phase 1)
    // -----------------------------------------------------------------------------------------------------------------
    inline CaseDef<2> make_case_2phase(const std::string& name, std::size_t axis, double mu_override)
    {
        using namespace materials;
        CaseDef<2> c;
        c.name          = name;
        Material copper = copper_2009(mu_override >= 0. ? mu_override : 9.2e10);
        c.mat           = {copper, air()};
        c.rho_init      = {8900., 1.};
        auto ax         = [axis](const Vec3& x) { return x[axis]; };
        auto tang       = [axis](double val)
        {
            Vec3 v{};
            v[(axis + 1) % 3] = val;
            return v;
        };
        auto norm = [axis](double val)
        {
            Vec3 v{};
            v[axis] = val;
            return v;
        };
        c.one_d = true;
        if (name == "advect" || name == "slip-advect") // Favrie 2009 sec. 5.1 and 6.1
        {
            c.eps   = 1e-7;
            c.tf    = 0.15e-3;
            const bool slip = name == "slip-advect";
            c.point = [=](const Vec3& x)
            {
                const bool s = ax(x) < 0.5;
                Vec3 v       = norm(1000.);
                if (slip && s)
                {
                    v[(axis + 1) % 3] = 100.;
                }
                return PointState{s ? 0u : 1u, v, 1e5};
            };
        }
        else if (name == "shear-slab") // Favrie 2009 sec. 6.2 (times corrected: 54.6 and 180 us)
        {
            c.eps   = 1e-7;
            c.tf    = 54.6e-6;
            c.times = {54.6e-6, 180e-6};
            c.point = [=](const Vec3& x)
            {
                const double s = ax(x);
                return PointState{(s > 0.15 && s < 0.85) ? 0u : 1u, tang((s > 0.35 && s < 0.65) ? 1000. : 0.), 1e5};
            };
        }
        else if (name == "double-expansion") // Favrie 2009 sec. 5.4
        {
            c.eps   = 1e-4;
            c.tf    = 88e-6;
            c.point = [=](const Vec3& x) { return PointState{0, norm(ax(x) < 0.5 ? -3000. : 3000.), 1e5}; };
        }
        else if (name == "solid-shear-tube") // Favrie 2009 sec. 5.5
        {
            c.eps   = 1e-6;
            c.tf    = 0.46e-4;
            c.point = [=](const Vec3& x)
            {
                const bool l = ax(x) < 0.5;
                return PointState{0, tang(l ? 0. : 100.), l ? 1e8 : 1e5};
            };
        }
        else if (name == "solid-gas-tube") // Favrie 2009 sec. 5.3 (t = 87 us assumed)
        {
            c.eps      = 1e-6;
            c.tf       = 87e-6;
            c.rho_init = {8900., 50.};
            c.point    = [=](const Vec3& x)
            {
                const bool l = ax(x) < 0.6;
                return PointState{l ? 0u : 1u, {}, l ? 5e9 : 1e5};
            };
        }
        else if (name == "impact1d") // copper slab at 800 m/s onto a copper slab at rest, in air
        {
            c.eps   = 1e-4;
            c.tf    = 40e-6;
            c.point = [=](const Vec3& x)
            {
                const double s = ax(x);
                return PointState{(s > 0.2 && s < 0.6) ? 0u : 1u, norm((s > 0.2 && s < 0.3) ? 800. : 0.), 1e5};
            };
        }
        else if (name == "wilkins") // Favrie & Gavrilyuk 2012 sec. 6.1: Wilkins' flying plate, aluminium, plastic
        {
            // mu printed as 248 GPa in the paper: 24.8 GPa assumed (aluminium)
            c.mat      = {Material{"aluminium", 3.4, 21.5e9, 24.8e9, 2700., 297.6e6}, air()};
            c.rho_init = {2700., 1.2};
            c.hi       = {0.05, 0.05, 0.05};
            c.base     = 0.0025;
            c.eps      = 1e-6;
            c.tf       = 5e-6;
            c.times    = {0.5e-6, 1e-6, 3e-6, 5e-6};
            c.plastic  = true;
            c.point    = [=](const Vec3& x)
            {
                const double s = ax(x);
                return PointState{s > 0.005 ? 0u : 1u, norm((s > 0.005 && s < 0.01) ? 800. : 0.), 1e5};
            };
        }
        else if (name == "impact2009") // Favrie 2009 sec. 7.1, elastic; --mu 9.2e10 (A), 1e9 (B), 0 (C)
        {
            c.one_d = false;
            c.eps   = 1e-4;
            c.level = 4;
            // output times of Fig. 16 (A, mu = 9.2e10), 17 (B, mu = 1e9), 18 (C, mu = 0)
            c.times = {36e-6, 106e-6, 320e-6, 609e-6};
            if (mu_override >= 0. && mu_override < 5e9)
            {
                c.times = mu_override > 0. ? std::vector<double>{35e-6, 140e-6, 420e-6, 710e-6}
                                           : std::vector<double>{75e-6, 187e-6, 600e-6, 1.03e-3};
            }
            c.tf = c.times.back();
            c.point = [](const Vec3& x)
            {
                const bool plate = x[0] > 0.2 && x[0] < 0.3 && x[1] > 0.25 && x[1] < 0.75;
                const bool proj  = x[0] > 0.1 && x[0] <= 0.2 && x[1] > 0.45 && x[1] < 0.55;
                return PointState{(plate || proj) ? 0u : 1u, Vec3{proj ? 800. : 0., 0., 0.}, 1e5};
            };
        }
        else if (name == "impact3d-disc" || name == "impact3d-square" || name == "impact3d-disc-plastic")
        {
            // 3D analogue of impact2009 (elastic) / impact2012 (plastic): a cylindrical projectile (axis x, diameter and
            // length 0.1 m, 800 m/s) against a circular (diameter 0.5 m) or square (0.5 m x 0.5 m) plate 0.1 m thick.
            // Domain 0.6 m cubed (use with hvi3d).
            c.one_d = false;
            c.lo    = {0., 0., 0.};
            c.hi    = {0.6, 0.6, 0.6};
            c.eps   = 1e-4;
            c.level = 2;
            c.tf    = 300e-6;
            const bool square  = name == "impact3d-square";
            const bool plastic = name == "impact3d-disc-plastic";
            const double gap   = plastic ? 0.005 : 0.;
            if (plastic)
            {
                c.mat     = {copper_2012(), air()};
                c.plastic = true;
            }
            c.point = [=](const Vec3& x)
            {
                const double r2   = sq(x[1] - 0.3) + sq(x[2] - 0.3); // axis of the projectile: y = z = 0.3
                const bool across = square ? (std::abs(x[1] - 0.3) < 0.25 && std::abs(x[2] - 0.3) < 0.25) : r2 < 0.25 * 0.25;
                const bool plate  = x[0] > 0.2 && x[0] < 0.3 && across;
                const bool proj   = x[0] > 0.1 - gap && x[0] <= 0.2 - gap && r2 < 0.05 * 0.05;
                return PointState{(plate || proj) ? 0u : 1u, Vec3{proj ? 800. : 0., 0., 0.}, 1e5};
            };
        }
        else if (name == "impact2012") // Favrie & Gavrilyuk 2012 sec. 6.2, elastic-plastic, narrow gap (positions assumed)
        {
            c.one_d    = false;
            c.mat      = {copper_2012(), air()};
            c.hi       = {0.7, 0.7, 0.7};
            c.eps      = 1e-4;
            c.level    = 4;
            c.plastic  = true;
            c.tf       = 1.1e-3;
            c.times    = {83e-6, 280e-6, 495e-6, 1.1e-3};
            c.point    = [](const Vec3& x)
            {
                const bool plate = x[0] > 0.3 && x[0] < 0.4 && x[1] > 0.1 && x[1] < 0.6;
                const bool proj  = x[0] > 0.195 && x[0] < 0.295 && x[1] > 0.3 && x[1] < 0.4;
                return PointState{(plate || proj) ? 0u : 1u, Vec3{proj ? 800. : 0., 0., 0.}, 1e5};
            };
        }
        else
        {
            throw std::invalid_argument("unknown case " + name);
        }
        return c;
    }

    // -----------------------------------------------------------------------------------------------------------------
    // cases with two solids (phases 0, 1) and air (phase 2)
    // -----------------------------------------------------------------------------------------------------------------
    inline CaseDef<3> make_case_3phase(const std::string& name, std::size_t axis = 0)
    {
        using namespace materials;
        CaseDef<3> c;
        c.name = name;
        if (name == "penetration2015") // Ndanou 2015 sec. 6.3: iron disc into an aluminium plate at 1000 m/s
        {
            c.mat      = {iron_2015(), aluminium_2015(), air()};
            c.rho_init = {7860., 2712., 1.};
            c.hi       = {0.08, 0.2, 0.2};
            c.base     = 0.005;
            c.eps      = 1e-6;
            c.level    = 4;
            c.plastic  = true;
            c.tf       = 100e-6;
            c.point    = [](const Vec3& x)
            {
                const bool disc  = sq(x[0] - 0.029) + sq(x[1] - 0.1) < 0.01 * 0.01;
                const bool plate = x[0] > 0.04 && x[0] < 0.06 && x[1] > 0.02 && x[1] < 0.18;
                return PointState{disc ? 0u : (plate ? 1u : 2u), Vec3{disc ? 1000. : 0., 0., 0.}, 1e5};
            };
        }
        else if (name == "kanel1d" || name == "spall1d") // Ndanou 2015 sec. 6.2.1 / 6.2.2: Al plate (2 mm, 700 m/s) onto Ti (9.8 mm)
        {
            // kanel1d: t = 0.2 us, analytical plateaus in Ndanou 2015 eq. (28)-(29); spall1d: t = 5 us (spall fractures)
            c.mat      = {aluminium_2015(), titanium_2015(), air()};
            c.rho_init = {2712., 4527., 1.};
            c.hi       = {0.017, 0.017, 0.017};
            c.base     = 0.001;
            c.eps      = 1e-6;
            c.level    = 7;
            c.plastic  = true;
            c.one_d    = true;
            c.tf       = name == "kanel1d" ? 0.2e-6 : 5e-6;
            c.point    = [axis](const Vec3& x)
            {
                const double s = x[axis];
                Vec3 v{};
                if (s > 0.001 && s < 0.003)
                {
                    v[axis] = 700.;
                }
                return PointState{(s > 0.001 && s < 0.003) ? 0u : ((s >= 0.003 && s < 0.0128) ? 1u : 2u), v, 1e5};
            };
        }
        else
        {
            throw std::invalid_argument("unknown case " + name);
        }
        return c;
    }

    inline bool is_3phase_case(const std::string& name)
    {
        return name == "penetration2015" || name == "kanel1d" || name == "spall1d";
    }
} // namespace hvi
