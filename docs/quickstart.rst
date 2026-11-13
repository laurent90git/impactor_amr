Quick start
===========

Build (conda environment with samurai 0.30, here ``cwisam_catapult``):

.. code-block:: bash

   conda activate cwisam_catapult
   cmake -S . -B build -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
   make -C build -j2            # hvi1d, hvi2d, hvi3d, test_physics
   ./build/test_physics         # checks of the physics kernel

Run:

.. code-block:: bash

   # 1D (uniform mesh, level 6 = 1280 cells on [0, 1])
   ./build/hvi1d --case shear-slab --level 6
   # 1D with multiresolution between levels 3 and 8
   ./build/hvi1d --case impact1d --level 8 --coarse-level 3
   # 2D plate impact of Favrie et al. (2009), case A, outputs at the times of Fig. 16
   ./build/hvi2d --case impact2009 --mu 9.2e10 --level 5 --coarse-level 2 --paper-times --tag A
   # 2D elastic-plastic plate impact of Favrie & Gavrilyuk (2012)
   ./build/hvi2d --case impact2012 --level 5 --coarse-level 2 --paper-times
   # 2D penetration of an iron disc into an aluminium plate, Ndanou et al. (2015)
   ./build/hvi2d --case penetration2015 --level 5 --coarse-level 2 --tf 120e-6 --nout 6

Main options (``--help`` for all):

=========================  ===========================================================================
``--case``                 test case (:doc:`verification`, :doc:`impacts`)
``--level``                finest level; cell size = base / 2^level (base: 0.05 m for most cases)
``--coarse-level``         coarsest level; below ``--level`` the mesh is adapted by multiresolution
``--adapt-eps``            multiresolution threshold (default 1e-3)
``--limiter``              ``minmod`` (default), ``vanleer``, ``mc``, ``none`` (first order)
``--rk``                   1 (forward Euler) or 2 (SSP-RK2, default)
``--cfl``                  CFL number on the sum of the wave speeds of all directions (default 0.4)
``--alim``                 fluidisation threshold of the slip treatment (default 1e-2, 0: none)
``--plastic``              1 / 0: plasticity on / off (default: the case's)
``--eps``                  residual volume fraction of the absent phases
``--mu``                   shear modulus of the solid of the two-phase cases
``--paper-times``          outputs at the times of the paper's figures
``--out``, ``--tag``       output directory and file suffix
=========================  ===========================================================================

Checkpoint and restart (samurai's HDF5 ``dump`` / ``load`` of the adapted mesh and the fields, plus the state of the
time loop):

.. code-block:: bash

   # a checkpoint every 15 min of wall time (default), and on SIGTERM / Ctrl-C (then the run stops cleanly)
   ./build/hvi2d --case impact2009 --level 5 --coarse-level 2 --paper-times --tag A --checkpoint-every 900
   kill -TERM <pid>            # -> out/impact2009_A_ckpt.h5 + .state
   # continue, with the same options plus --resume; outputs and log continue where they stopped
   ./build/hvi2d --case impact2009 --level 5 --coarse-level 2 --paper-times --tag A --resume out/impact2009_A_ckpt

A resumed run is bitwise identical to an uninterrupted one (checked in 1D and 2D, adapted meshes).

Outputs, in ``<out>/<case>[_<tag>]``:

- ``_NNN.npy``: one row per cell, columns listed in ``_columns.txt`` (centre, size, level, density, solid fraction,
  velocity, pressure, stress components, ``A[1][0]`` of the first solid, slip switch, von Mises stress);
- ``.log``: per logged step, time, number of cells, masses of the phases, momentum, total energy, counters of the
  relaxation.

Plots: ``scripts/schlieren.py`` (2D: schlieren, or any field with ``--field``), ``scripts/compare_proto.py`` (1D:
comparison with the NumPy prototype ``prototype/hvi1d.py``).
