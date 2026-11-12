// Material parameters and the hydrodynamic (volumetric) equation of state: stiffened gas.
//   rho e^h = (p + gamma p_inf) / (gamma - 1),   c^2 = gamma (p + p_inf) / rho
// A material with mu = 0 is a fluid. yield <= 0 means purely elastic.
#pragma once

#include <cmath>
#include <string>

namespace hvi
{
    struct Material
    {
        std::string name;
        double gamma = 1.4;  ///< stiffened-gas exponent
        double pinf  = 0.;   ///< stiffened-gas pressure [Pa]
        double mu    = 0.;   ///< shear modulus [Pa] (0: fluid)
        double rho0  = 1.;   ///< reference density of the shear energy [kg/m3]
        double yield = 0.;   ///< von Mises yield stress [Pa] (<= 0: elastic)

        [[nodiscard]] bool is_solid() const
        {
            return mu > 0.;
        }

        [[nodiscard]] bool is_plastic() const
        {
            return mu > 0. && yield > 0.;
        }

        /// volumetric hydrodynamic energy rho e^h of the phase (per unit volume of the phase)
        [[nodiscard]] double rho_eh(double p) const
        {
            return (p + gamma * pinf) / (gamma - 1.);
        }

        /// pressure from rho e^h
        [[nodiscard]] double pressure(double rho_eh_) const
        {
            return (gamma - 1.) * rho_eh_ - gamma * pinf;
        }

        /// squared hydrodynamic sound speed
        [[nodiscard]] double c2(double rho, double p) const
        {
            return gamma * (p + pinf) / rho;
        }

        /// lowest admissible pressure (p > -p_inf)
        [[nodiscard]] double pmin() const
        {
            return -pinf;
        }
    };

    // Parameters used by the reference papers -------------------------------------------------------------------------
    namespace materials
    {
        /// Favrie, Gavrilyuk, Saurel (2009): copper with the shear modulus of copper at ambient temperature
        inline Material copper_2009(double mu = 9.2e10)
        {
            return {"copper", 4.22, 342.0e8, mu, 8900., 0.};
        }

        /// Favrie & Gavrilyuk (2012), section 6.2: "copper" alloy, elastic-plastic
        inline Material copper_2012()
        {
            return {"copper", 4.22, 34.2e9, 43.3e9, 8900., 600e6};
        }

        /// Ndanou, Favrie, Gavrilyuk (2015), table 1
        inline Material iron_2015()
        {
            return {"iron", 3.9, 43.6e9, 82e9, 7860., 200e6};
        }

        inline Material aluminium_2015()
        {
            return {"aluminium", 3.5, 32e9, 26e9, 2712., 60e6};
        }

        inline Material titanium_2015()
        {
            return {"titanium", 2.6, 44e9, 42e9, 4527., 1030e6};
        }

        inline Material air()
        {
            return {"air", 1.4, 0., 0., 1., 0.};
        }
    } // namespace materials
} // namespace hvi
