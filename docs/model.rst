Physical model
==============

The model is the multi-solid, multi-fluid diffuse-interface model of Ndanou, Favrie & Gavrilyuk [NFG15]_, which
extends the solid–fluid model of Favrie, Gavrilyuk & Saurel [FGS09]_ and Favrie & Gavrilyuk [FG12]_. All the
materials share one velocity; each phase has its own volume fraction, density, entropy (internal energy) and, for
solids, its own elastic deformation. Interfaces are regions where several volume fractions are non-zero; a residual
fraction :math:`\varepsilon` of every phase is present everywhere (1e-4 to 1e-7).

Kinematics of a solid phase
---------------------------

Let :math:`X^\beta` be the Lagrangian coordinates of solid :math:`l` and :math:`e^\beta = \nabla_x X^\beta` its
cobasis vectors. The code stores the matrix :math:`A = F^{-1}` whose row :math:`\beta` is :math:`e^\beta`
(``A[b][i] = dX^b/dx_i``). The Finger tensor and its unimodular part are

.. math::

   G = A^T A = F^{-T} F^{-1}, \qquad g = \frac{G}{|G|^{1/3}} = \frac{G}{|\det A|^{2/3}}.

Without plasticity the cobasis is transported by

.. math::
   :label: cobasis

   \frac{D e^\beta}{Dt} + \left(\frac{\partial v}{\partial x}\right)^T e^\beta = 0.

With plasticity, :eq:`cobasis` gets a relaxation source (below) and describes the *elastic* deformation.

Energy and stress
-----------------

The specific energy of each phase is in separable form [GFS08]_ [NFG15]_: a hydrodynamic part depending on density
and entropy, and a shear part unaffected by volume changes,

.. math::

   e_l = e^h_l(\rho_l, \eta_l) + e^e_l(g_l), \qquad
   \rho_l e^h_l = \frac{p_l + \gamma_l p_{\infty,l}}{\gamma_l - 1}, \qquad
   e^e_l = \frac{\mu_l}{8 \rho_{0,l}} \left(\operatorname{tr}(g_l^2) - 3\right).

The pressure comes only from the stiffened-gas part, the deviatoric stress only from the shear part:

.. math::

   \sigma_l = -p_l I + S_l, \qquad
   S_l = -2\rho_l \frac{\partial e^e_l}{\partial G_l} G_l
       = -\frac{\mu_l \rho_l}{2 \rho_{0,l}}\left(g_l^2 - \frac{\operatorname{tr}(g_l^2)}{3} I\right).

A fluid is a phase with :math:`\mu_l = 0`. For small strains Hooke's law with shear modulus :math:`\mu_l` is
recovered. The mixture stress is :math:`\sigma = -\sum_l \alpha_l p_l\, I + \sum_l \alpha_l S_l`.

Governing equations
-------------------

The pressure non-equilibrium model (Ndanou 2015, eq. 15–16) solved by the code is

.. math::

   \frac{D\alpha_l}{Dt} &= \lambda (p_l - p_I),\\
   \frac{\partial \alpha_l\rho_l}{\partial t} + \nabla\cdot(\alpha_l \rho_l v) &= 0,\\
   \frac{\partial \rho v}{\partial t} + \nabla\cdot(\rho v\otimes v - \sigma) &= 0,\\
   \frac{\partial \rho E}{\partial t} + \nabla\cdot(\rho E v - \sigma v) &= 0,\\
   \frac{\partial \alpha_l\rho_l e^h_l}{\partial t} + \nabla\cdot(\alpha_l\rho_l e^h_l v) + \alpha_l p_l \nabla\cdot v &= -p_I \frac{D\alpha_l}{Dt},\\
   \frac{D e^\beta_l}{Dt} + \left(\frac{\partial v}{\partial x}\right)^T e^\beta_l &= \frac{1}{\tau_l} S_l e^\beta_l,

with :math:`\rho E = \sum_l \alpha_l\rho_l (e^h_l + e^e_l) + \rho |v|^2/2`. The pressure relaxation rate
:math:`\lambda` is infinite: the pressures are equalised after each step (the interfacial pressure :math:`p_I` then
drops out). The plastic relaxation time :math:`\tau_l` is zero above the yield surface (ideal plasticity) and
infinite below it.

The state vector of a cell, for N phases of which the first NS are solids, is

=========================  ==========  =====================================
:math:`\alpha_l`           N           non-conservative
:math:`\alpha_l \rho_l`    N           conservative
:math:`\rho v`             3           conservative (always 3 components)
:math:`\rho E`             1           conservative
:math:`\alpha_l\rho_l e_l` N           non-conservative
:math:`A_s`                9 NS        non-conservative
=========================  ==========  =====================================

so 19 components for one solid and air, 31 for two solids and air (``physics/model.hpp``, ``Layout``).

Wave speeds
-----------

For a direction :math:`d` the acoustic tensor (per unit mass) of the mixture is

.. math::

   K = c_f^2\, e_d e_d^T + \sum_{s\ \mathrm{solid}} Y_s D_s, \qquad c_f^2 = \sum_l Y_l \gamma_l\frac{p_l + p_{\infty,l}}{\rho_l},

.. math::

   D_s = \frac{\mu_s}{2\rho_{0,s}}\Big[g_{dd}\, g + (g e_d)(g e_d)^T + g^2
         - \tfrac13\big((g^2 e_d) e_d^T + e_d (g^2 e_d)^T\big) + \tfrac{\operatorname{tr} g^2}{9} e_d e_d^T\Big],

derived for this code from the second derivatives of the shear energy (checked against the uniaxial value
:math:`\frac{\mu}{\rho_0}(\frac{11}{9} a^{8/3} + \frac19 a^{-4/3})` for :math:`A = \mathrm{diag}(a,1,1)` in
``tests/test_physics.cpp``). At rest it gives :math:`4\mu/3\rho_0` along :math:`d` and :math:`\mu/\rho_0` across.
The longitudinal speed used by the Riemann solver is the Gershgorin bound of the largest eigenvalue of :math:`K`, the
transverse one :math:`c_t^2 = \max(K_{t_1t_1}, K_{t_2t_2}) + |K_{t_1t_2}|`.

Plasticity
----------

Von Mises criterion :math:`S:S \le \frac23 Y^2`. Following Godunov & Romenskii, Favrie & Gavrilyuk [FG12]_ §3.4 and
Ndanou et al. §5.2, the relaxation acts on the singular values of the cobasis with frozen singular vectors. With
:math:`G = U \operatorname{diag}(\lambda_a) U^T`, :math:`b_a = \lambda_a / (\lambda_1\lambda_2\lambda_3)^{1/3}` and
:math:`x_a = \ln b_a`, the principal deviatoric stresses are
:math:`S_a = -\frac{\mu\rho}{2\rho_0}\left(b_a^2 - \frac13\sum_c b_c^2\right)`, and the relaxation path is

.. math::

   \frac{dx_a}{ds} = \frac{S_a}{|S|}, \qquad \sum_a x_a = 0 \ \text{(density unchanged)},

followed until :math:`S:S = \frac23 Y^2` (``physics/plasticity.hpp``). The cobasis becomes
:math:`A \leftarrow A\, U \operatorname{diag}(e^{\Delta x_a / 2}) U^T`. The released elastic energy turns into heat
through the energy reset of the relaxation step: total energy is conserved and entropy increases.

Cavitation and fracture
-----------------------

Fracture is not a separate model. The residual gas in a solid has a non-negative pressure, so the relaxed mixture
pressure cannot go below zero: in strong tension the gas fraction grows (void opening), which is how spall and
cracks appear [NFG15]_ §6. No damage parameter is involved; results in tension therefore depend on the residual
fraction and on the mesh.

Slip at solid/fluid interfaces
------------------------------

With residual fractions, the "fluid" side of an interface contains some solid, which has a shear stiffness: without
treatment the solid drags the fluid (no slip) and sends shear waves into it. Favrie et al. correct this with a
level-set and a ghost-fluid treatment of the tangential velocity [FGS09]_ §5.6–6, [FG12]_ §5.2. This code uses a
local rule instead (``physics/slip.hpp``, policy ``Fluidisation``):

.. math::

   \chi = \begin{cases} 1 & \sum_s \alpha_s \ge \alpha_{lim}\\ 0 & \text{otherwise}\end{cases}

Where :math:`\chi = 0` the solids have no shear stress, no shear energy and no shear impedance, and their cobasis is
reset to a pure dilatation :math:`A \leftarrow |\det A|^{1/3} I`. Combined with the shear Riemann solver, which
weights the tangential traction by the shear impedances, a face between a solid cell and a fluid cell carries zero
tangential traction: free slip with a single conservative flux per face. Default :math:`\alpha_{lim} = 10^{-2}`;
see :doc:`verification` for the effect of the threshold.
