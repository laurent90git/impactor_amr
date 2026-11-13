Numerical method
================

One time step
-------------

The equations are split into a hyperbolic part :math:`L` (fluxes and non-conservative products, all directions
together) and a cell-local relaxation :math:`R`. Explicit SSP-RK2:

.. math::

   U^{(1)} = R\big(U^n - \Delta t\, L(U^n)\big), \qquad
   U^{n+1} = R\Big(\tfrac12\big(U^n + U^{(1)} - \Delta t\, L(U^{(1)})\big)\Big),

with :math:`\Delta t = \mathrm{CFL}\; h_{min} / \max_{cells} \sum_d (|v_d| + c_{l,d})`, CFL = 0.4. ``--rk 1`` gives
forward Euler. The papers use MUSCL-Hancock; SSP-RK2 maps directly onto samurai's scheme operators.

When the mesh is adaptive, it is adapted before each step (multiresolution, below) and :math:`R` is applied to the
new cells.

Space operator
--------------

For each face :math:`j+1/2` of direction :math:`d`, from the four cells :math:`j-1 \dots j+2`
(``solver/flux_operator.hpp``):

1. conversion of the four cells to primitive variables
   :math:`(\alpha_l, \rho_l, v, p_l, A_s)` (``Model::to_prim``);
2. MUSCL reconstruction with a limiter (minmod by default) of all of them; the slip switch :math:`\chi` is not
   reconstructed, each face state keeps that of its cell; a non-admissible face state falls back to the cell value
   (``numerics/reconstruction.hpp``);
3. face solver (``numerics/riemann.hpp``) returning the flux :math:`F` and the face velocity :math:`v^*`;
4. the contributions to the two cells, including the non-conservative products evaluated with each cell's own values.

The cell update is

.. math::

   \frac{dU_i}{dt} = -\frac{F_{i+1/2} - F_{i-1/2}}{h} - C(U_i)\,\frac{v^*_{i+1/2} - v^*_{i-1/2}}{h},

with the non-conservative products (one direction :math:`d`, :math:`t \ne d` the tangential directions)

.. math::

   \partial_t \alpha_l &\mathrel{+}= \alpha_l\, \partial_d u^*_d, \qquad
   \partial_t (\alpha\rho e)_l \mathrel{-}= \alpha_l p_l\, \partial_d u^*_d,\\
   \partial_t A_{b,t} &\mathrel{+}= A_{b,t}\, \partial_d u^*_d, \qquad
   \partial_t A_{b,d} \mathrel{-}= \textstyle\sum_{t\ne d} A_{b,t}\, \partial_d v^*_t .

Summed over the directions these give the multi-dimensional :eq:`cobasis`. In samurai this is a non-linear
flux-based scheme with a *non-conservative* flux function: the contribution to the left cell of a face is
:math:`F + NC(U_j, v^*)`, to the right cell :math:`-(F + NC(U_{j+1}, v^*))`.

Face solver
-----------

``LongitudinalHllcShear`` sums two solvers, after the splitting idea of Favrie, Gavrilyuk & Ndanou [FGN14]_ (which
applies the sub-systems one after the other; here they are summed in one flux):

*Longitudinal HLLC* (mass, normal momentum, energy; tangential momentum only advected). Davis wave speeds with the
longitudinal mixture speed, contact speed

.. math::

   S_M = \frac{\sigma_{n,L} - \sigma_{n,R} + \rho_L u_L(S_L - u_L) - \rho_R u_R (S_R - u_R)}{\rho_L(S_L - u_L) - \rho_R(S_R - u_R)},

star states of the upwind side :math:`K` with the factor :math:`(S_K - u_K)/(S_K - S_M)` on the densities and on
the normal column of the cobasis; phase energies along their Hugoniot (Saurel et al. [SPB09]_ eq. II.8), else their
isentrope.

*Shear solver* for each tangential direction :math:`t`, with shear impedances :math:`Z = \rho c_t`:

.. math::

   v^*_t = \frac{Z_L v_{t,L} + Z_R v_{t,R} + \tau_R - \tau_L}{Z_L + Z_R}, \qquad
   \tau^* = \frac{Z_R\tau_L + Z_L\tau_R + Z_L Z_R (v_{t,R} - v_{t,L})}{Z_L + Z_R},

added to the tangential momentum and energy fluxes. With :math:`Z = 0` on one side the traction vanishes: this is
what makes the slip treatment work with one flux per face.

Relaxation step
---------------

``physics/relaxation.hpp``, cell by cell:

1. volume fractions clipped to :math:`[10^{-12}, 1]` and normalised;
   **vacuum guard**: the scheme is not positivity-preserving, and in the strongest expansions (corner jets) a cell
   can be drained to vacuum. If the mixture density falls below :math:`\rho_{floor}` (1e-3 times the lightest
   initial density) the lightest phase is added up to :math:`\rho_{floor}`, the velocity is capped at 20 km/s and
   the pressure set to the floor (counted ``vacuum``; mass and energy are not conserved in these cells). The same
   reset is applied when a stiff phase is drained below 1e-3 of its reference density while keeping its volume
   fraction (its volume goes to the gas): its stiffened-gas sound speed would otherwise diverge;
2. slip policy: reset of the cobasis of each solid whose volume fraction is below :math:`\alpha_{lim}`;
   **failure safeguard**: the same reset where the elastic strain is extreme, :math:`\operatorname{tr} g > 50`
   (principal stretch about 7), which only happens in smeared interface cells of jets (counted ``failed``);
3. plastic relaxation of each solid to the yield surface;
4. instantaneous pressure relaxation (Saurel et al. [SPB09]_ §3.3): solve
   :math:`\sum_l \alpha_l \frac{p_{0,l} + \gamma_l p_{\infty,l} + (\gamma_l - 1)p}{\gamma_l (p + p_{\infty,l})} = 1`
   by Newton in :math:`p > \max_l(-p_{\infty,l})`, update the volume fractions;
5. energy reset from the total energy with the mixture stiffened-gas law (Ndanou 2015 eq. 26). In strong expansion
   the stiff phases (copper) can be drained while keeping their volume fraction: the hydrodynamic energy is then too
   low for any admissible pressure. **Void opening**: at the floor pressure the volume fractions of the stiff phases
   are scaled by :math:`k < 1` and the gaseous phases fill the rest, :math:`k` being chosen so that the mixture
   equation of state matches the total energy exactly (counted ``voided``). Only when this is impossible is the
   relaxed pressure kept and the total energy changed (counted ``floored``).

The first version of the code (cases A and C of :doc:`impacts`) kept the relaxed pressure in step 5, i.e. it injected
energy; in case B this produced hot, very light copper with a sound speed of 65 km/s and the time step collapsed.
The diagnostics printed by the solver when the time step collapses (state of the limiting cell) found the three
problems fixed by the guards above.

Multiresolution
---------------

samurai's multiresolution adaptation (``mra_config().epsilon(1e-3).regularity(1)``) on an indicator with the total
solid fraction, :math:`\ln\rho / \ln(\rho_{max}/\rho_{min})`, the pressure and the velocity normalised by reference
values of the case. Interfaces and shocks stay at the finest level. Ghosts at level jumps are predicted by samurai;
the face solver works on them as on any state (``to_prim`` clips the volume fractions).

Boundary conditions
-------------------

Transmissive on all sides: zero-order extrapolation into the two ghost layers (``solver/boundary.hpp``).
