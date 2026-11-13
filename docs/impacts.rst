2D impacts
==========

Figures produced by ``scripts/make_doc_figures_2d.py`` from the runs of ``scripts/run_production.sh``: second order,
multiresolution between levels 2 and 5, transmissive boundaries. The schlieren images show
:math:`\exp(-k |\nabla \ln\rho| / \max|\nabla \ln\rho|)` on the finest grid.

.. include:: _static/fig/impacts_table.rst

The mass changes only through the outflow boundaries. The energy changes come from the outflow and from the cells
where the relaxation could not keep an admissible pressure (``floored``, see :doc:`numerics`).

Copper projectile on a copper plate, elastic (Favrie et al. 2009, §7.1)
-----------------------------------------------------------------------

Domain 1 m x 1 m, plate [0.2, 0.3] x [0.25, 0.75] m, projectile [0.1, 0.2] x [0.45, 0.55] m at 800 m/s, in contact
at t = 0. Copper :math:`\gamma = 4.22`, :math:`p_\infty = 34.2` GPa, :math:`\rho = 8900` kg/m³; air :math:`\gamma =
1.4`, 1 kg/m³; residual volume fraction 1e-4. Finest cell 1.56 mm (640 x 640 equivalent; the paper uses
1000 x 1000).

Case A, :math:`\mu = 92` GPa (Fig. 16 of the paper, t = 36, 106, 320, 609 µs):

.. figure:: _static/fig/impacts_A_placeholder.png
   :width: 100%

Case B, :math:`\mu = 1` GPa (Fig. 17, t = 35, 140, 420, 710 µs):

.. figure:: _static/fig/impacts_B_placeholder.png
   :width: 100%

Case C, :math:`\mu = 0` (fluid limit, Fig. 18, t = 75, 187, 600, 1030 µs):

.. figure:: _static/fig/impacts_C_placeholder.png
   :width: 100%

Elastic-plastic copper (Favrie & Gavrilyuk 2012, §6.2)
-------------------------------------------------------

Domain 0.7 m x 0.7 m, copper with :math:`\mu = 43.3` GPa, :math:`Y = 600` MPa, projectile at 800 m/s separated
from the plate by a 5 mm gap (positions assumed, see :doc:`references`). Paper times 83, 280, 495 µs and 1.1 ms
(Fig. 7: the plate breaks into five pieces).

.. figure:: _static/fig/impacts_P_placeholder.png
   :width: 100%

Iron disc into an aluminium plate (Ndanou et al. 2015, §6.3)
------------------------------------------------------------

Three materials: iron disc (radius 1 cm, 1000 m/s), aluminium plate 2 cm x 16 cm, air; both metals elastic-plastic
(table 1 of the paper). Domain 8 cm x 20 cm (assumed), finest cell 0.16 mm.

.. figure:: _static/fig/impacts_N_placeholder.png
   :width: 100%
