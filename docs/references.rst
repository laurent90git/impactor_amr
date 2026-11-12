References
==========

.. [FGS09] N. Favrie, S. L. Gavrilyuk, R. Saurel. *Solid–fluid diffuse interface model in cases of extreme
   deformations.* Journal of Computational Physics 228 (2009) 6037–6077. doi:10.1016/j.jcp.2009.05.015

.. [FG12] N. Favrie, S. L. Gavrilyuk. *Diffuse interface model for compressible fluid – compressible elastic–plastic
   solid interaction.* Journal of Computational Physics 231 (2012) 2695–2723. doi:10.1016/j.jcp.2011.11.027

.. [FGN14] N. Favrie, S. Gavrilyuk, S. Ndanou. *A thermodynamically compatible splitting procedure in
   hyperelasticity.* Journal of Computational Physics 270 (2014) 300–324. doi:10.1016/j.jcp.2014.03.051

.. [NFG15] S. Ndanou, N. Favrie, S. Gavrilyuk. *Multi-solid and multi-fluid diffuse interface model: applications to
   dynamic fracture and fragmentation.* Journal of Computational Physics 295 (2015) 523–555.
   doi:10.1016/j.jcp.2015.04.024

.. [SPB09] R. Saurel, F. Petitpas, R. A. Berry. *Simple and efficient relaxation methods for interfaces separating
   compressible fluids, cavitating flows and shocks in multiphase mixtures.* Journal of Computational Physics 228
   (2009) 1678–1712.

.. [GFS08] S. L. Gavrilyuk, N. Favrie, R. Saurel. *Modelling wave dynamics of compressible elastic materials.*
   Journal of Computational Physics 227 (2008) 2941–2969.

.. [WBN21] T. Wallis, P. T. Barton, N. Nikiforakis. *A flux-enriched Godunov method for multi-material problems with
   interface slide and void opening.* Journal of Computational Physics (2021), arXiv:2011.05569. (An alternative
   treatment of slip in diffuse-interface methods, not implemented.)

samurai: adaptive multiresolution library, https://github.com/hpc-maths/samurai (version 0.30.1 used here).

Notes on the papers
-------------------

- Favrie 2009 §6.2: the times printed in the text and in Fig. 14 (5.46 and 18 µs) are a factor 10 too small; the
  wave fronts correspond to 54.6 and 180 µs. The axes of Fig. 14 are mutually inconsistent (a 500 m/s plateau
  implies a shear strain of 0.155 and a shear stress of 14.3 GPa).
- Favrie 2009 §5.3: the time of 870 µs is too long for the 1 m domain; 87 µs is used here.
- Favrie 2009 §7.1: positions of the plate and of the projectile measured on Fig. 16 (the sketch of Fig. 15 is not
  to scale): plate [0.2, 0.3] x [0.25, 0.75] m, projectile [0.1, 0.2] x [0.45, 0.55] m.
- Favrie 2012 §6.1: the shear modulus of aluminium is printed as 248 GPa; 24.8 GPa is used.
- Favrie 2012 §6.2: the positions and the width of the gap are not given; plate [0.3, 0.4] x [0.1, 0.6] m,
  projectile [0.195, 0.295] x [0.3, 0.4] m (gap 5 mm) are used.
- Ndanou 2015 §6.3: the domain size is not given; [0, 0.08] x [0, 0.2] m, plate [0.04, 0.06] x [0.02, 0.18] m,
  disc centred at (0.029, 0.1) m (gap 1 mm) are used.
