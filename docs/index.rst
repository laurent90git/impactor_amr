Eulerian hypervelocity impact solver
====================================

An Eulerian solver for impacts of elastic(-plastic) solids in a compressible fluid, with diffuse interfaces, built on
the adaptive multiresolution library `samurai <https://github.com/hpc-maths/samurai>`_. It implements the
multi-material hyperelastic model of Ndanou, Favrie & Gavrilyuk (2015) with an explicit Godunov-type scheme and
reproduces the impact tests of Favrie, Gavrilyuk & Saurel (2009), Favrie & Gavrilyuk (2012) and Ndanou et al. (2015).

.. toctree::
   :maxdepth: 2

   quickstart
   model
   numerics
   code
   verification
   impacts
   references

Status
------

- Physics: N phases (solids and fluids), hyperelasticity in separable form, ideal von Mises plasticity, cavitation
  and fracture through void opening, slip at solid/fluid interfaces.
- Numerics: MUSCL + SSP-RK2, longitudinal HLLC + shear impedance solver, instantaneous pressure and plastic
  relaxations, multiresolution adaptation (samurai).
- Verified in 1D against an independent NumPy implementation (round-off agreement) and against analytical plateau
  values of Ndanou et al. (2015); 2D impact tests in :doc:`impacts`.
