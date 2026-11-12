Code organisation
=================

::

   src/
     core/linalg.hpp              3x3 algebra, symmetric eigen-decomposition (Jacobi)
     physics/                     no samurai inside: testable, reusable
       material.hpp               Material (stiffened gas + shear modulus + yield), parameters of the papers
       shear_energy.hpp           NeoHookeanSeparable: energy, stress, acoustic tensor
       slip.hpp                   Fluidisation, NoSlipTreatment
       model.hpp                  Layout<N,NS>, Prim, Mixture, Model: conversions, mixture stress, wave speeds
       plasticity.hpp             VonMisesInstantaneous, NoPlasticity
       relaxation.hpp             CellRelaxation: slip reset, plasticity, pressure relaxation, energy reset
     numerics/
       reconstruction.hpp         MusclPrimitive and the limiters
       riemann.hpp                LongitudinalHllcShear (face solver)
     solver/                      samurai glue
       flux_operator.hpp          make_flux_operator: reconstruction + face solver as a samurai scheme
       boundary.hpp               transmissive boundary condition
       run.hpp                    run_case: mesh, initial condition, time loop, adaptation, outputs
     cases/cases.hpp              test cases
     io/npy.hpp                   .npy writer
     app/main.cpp                 command line
   tests/test_physics.cpp         kernel checks
   prototype/hvi1d.py             independent 1D NumPy implementation (reference)
   scripts/                       plots, comparison with the prototype, production runs

Flow of one run
---------------

``main.cpp`` parses the options and picks the model instantiation from the case (``Model<2,1>``: one solid and air;
``Model<3,2>``: two solids and air), then calls ``run_case<dim>``:

1. build the samurai mesh (``mra::make_mesh``) and the fields;
2. ``make_flux_operator`` builds :math:`L` from the model, a reconstruction and a face solver;
3. initial condition by sub-sampling the case's point description (``initialise``), mesh adapted to it;
4. time loop: adapt, time step from the wave speeds, SSP-RK2 stages each followed by ``relax_all``, diagnostics,
   outputs.

Extension points
----------------

Every building block is a small struct with one documented interface, passed as a template parameter:

=======================  ==========================================================  =============================
Block                    Interface                                                   Default
=======================  ==========================================================  =============================
shear energy             ``evaluate(A)``, ``energy``, ``deviator``, ``acoustic``,     ``NeoHookeanSeparable``
                         ``principal_deviator``
slip treatment           ``chi(alpha_solid)``, ``reset(A)``                           ``Fluidisation``
plasticity               ``bool operator()(Mat3& A, const Material&, double rho)``    ``VonMisesInstantaneous``
reconstruction           ``operator()(model, 4 primitive states, wl, wr)``            ``MusclPrimitive``
face solver              ``FaceResult operator()(model, wl, wr, d)``                  ``LongitudinalHllcShear``
=======================  ==========================================================  =============================

To try another Riemann solver, write a struct with the face-solver interface returning the flux of every component
(the advective part only for the non-conservative ones) and the face velocity, and pass it to
``make_flux_operator<Field>(model, recon, solver)``. The non-conservative products are added by the operator, the
relaxation is unchanged. To change the time integration, edit the loop in ``run.hpp`` (the operator ``op(U)``
returns the divergence; ``relax_all`` is the cell-local step).

Testing a change
----------------

1. ``./build/test_physics``;
2. 1D ladder against the prototype: ``./build/hvi1d --case <case> --level 4 --out out/t`` then
   ``python scripts/compare_proto.py <case> out/t/<case>_000.npy --n 320`` (differences must stay at round-off as long
   as the prototype implements the same method; otherwise update the prototype first);
3. a 2D run with the conservation columns of the log.
