# Development plan: Eulerian hyperelastic diffuse-interface solver on Samurai

Audience: the developer (human or model) who implements the solver. Read this file fully before writing code.
Everything stated as "validated" was run in the 1D prototype `prototype/hvi1d.py`; everything else is design and is
marked as such.

## 1. Goal and definition of done

Replicate section 7.1 of Favrie, Gavrilyuk, Saurel, JCP 228 (2009) ("Impact of a projectile on a solid plate", Fig. 15-18)
with a C++ solver built on Samurai (multiresolution adaptation, explicit time integration), using the model of
Ndanou, Favrie, Gavrilyuk, JCP 295 (2015).

Done means:

1. The 1D verification ladder of section 7 passes on uniform and adapted meshes.
2. The 2D impact runs to the last time of the paper for the three shear moduli and produces density-schlieren images
   at the paper's times that show the same qualitative behaviour (section 2).
3. The adapted run agrees with a uniform run at the same finest level (interface position, plate shape) and uses
   clearly fewer cells.
4. A `README.md` gives build and run commands and the results obtained.

Out of scope for the first delivery, but the code must not block them: plasticity, more than one solid, 3D, MPI.

## 2. Target test (Favrie 2009, section 7.1)

| Item | Value |
|---|---|
| Domain | 1 m x 1 m; paper mesh 1000 x 1000 uniform |
| Plate (at rest) | 0.1 m (x) by 0.5 m (y): `[0.2, 0.3] x [0.25, 0.75]` |
| Projectile | 0.1 m square, `[0.1, 0.2] x [0.45, 0.55]`, velocity (800, 0) m/s, already in contact at t = 0 |
| Copper | stiffened gas gamma = 4.22, p_inf = 342e8 Pa, rho = 8900 kg/m3 |
| Shear modulus | case A: 9.2e10 Pa; case B: 1.0e9 Pa; case C: 0 |
| Air | ideal gas gamma = 1.4, rho = 1 kg/m3, p = 1e5 Pa everywhere at t = 0 |
| Residual fractions | 1e-4 of air in copper, 1e-4 of copper in air |
| Scheme in the paper | MUSCL-Hancock, minmod, HLLC |
| Output | schlieren of the mixture density |
| Times A | 36, 106, 320, 609 microseconds |
| Times B | 35, 140, 420, 710 microseconds |
| Times C | 75, 187, 600, 1030 microseconds |

Notes:

- The positions come from measuring Fig. 16 (first panel, where the initial plate outline is drawn). The sketch of
  Fig. 15 is not to scale. Sizes are from the text.
- Boundary conditions are not given in the paper. Use transmissive conditions on all sides.
- Projectile and plate are the same material and the same phase, so they are welded on contact.
- Expected behaviour: A, the T-shaped body stays in one piece and the plate oscillates; B, the plate bends strongly
  but the strain stays bounded; C, fluid limit, the plate keeps deforming and filaments form.
- Check value: momentum conservation gives a final mean velocity of 800/6 = 133 m/s for the body.
- With 1e-4 of copper in the air, the "air" is 47 % copper by mass (mixture density 1.89 kg/m3). This is the paper's
  choice; keep it for the replication and make it a parameter.

Follow-up target once this works: Favrie & Gavrilyuk, JCP 231 (2012), section 6.2. Same geometry with plasticity
(domain 0.7 m, 700 x 700, mu = 43.3 GPa, yield 600 MPa, a narrow gap between projectile and plate, times 83, 280, 495
microseconds and 1.1 ms). The plate breaks into five pieces.

## 3. Reference material

| File | Use |
|---|---|
| `Ndanou et al. - 2015 ...pdf` | the model (sections 2, 4, 5) |
| `Favrie et al. - 2009 ...pdf` | the target test (7.1), 1D tests (5.1 to 6.2), MUSCL details (appendix B) |
| `favrie2012.pdf` | slip treatment of the authors (5.2), plasticity (3.4, 5.3), follow-up test (6.2) |
| `favrie2014.pdf` | splitting into a longitudinal and shear sub-systems (sections 3, 4), basis of our flux |
| `Richard Saurel et al. - 2009 ...pdf` | pressure-relaxation method for the fluid case |
| `prototype/hvi1d.py` | executable reference of every formula of sections 4 and 5, and of the 1D tests |
| `~/GIT/CREATION_LIBRAIRIE_COUPLAGE/new_library_dev/examples/samurai_euler/` | Samurai patterns to copy: `euler2d.cpp` (mesh, fields, time loop, npy output), `weno_euler.hpp` (flux-based scheme), `bc_euler.hpp` (boundary conditions filling R ghost layers), `turek_fsi_fluid.cpp` lines 705-730 (multiresolution adaptation), `CMakeLists.txt` |

Known errors in Favrie 2009 that affect the 1D tests:

- Section 6.2: the times are 54.6 and 180 microseconds, not 5.46 and 18. The moving region is `[0.35, 0.65]`.
  The axes of Fig. 14 are mutually inconsistent (a 500 m/s plateau implies a shear strain of 0.155 and 14.3 GPa).
- Section 5.3: the quoted time of 870 microseconds is too long for a 1 m domain; 87 microseconds is assumed.

## 4. Model

### 4.1 State

N phases share one velocity. The state is always three-dimensional, whatever the mesh dimension: three velocity
components and a 3 x 3 matrix per solid phase. This lets the 1D tests with transverse velocity, the 2D case and 3D
use the same kernel.

Per cell, for phases l = 1..N and solid phases s = 1..Ns (phases with mu > 0):

| Symbol | Count | Meaning | Type of equation |
|---|---|---|---|
| alpha_l | N | volume fraction | non-conservative |
| (alpha rho)_l | N | partial density | conservative |
| rho v | 3 | momentum | conservative |
| rho E | 1 | total energy | conservative |
| (alpha rho e^h)_l | N | hydrodynamic internal energy | non-conservative |
| A_s | 9 Ns | A = F^-1, `A[b][i] = dX^b/dx_i` (row b is the cobasis vector e^b) | non-conservative |

For copper and air: 2 + 2 + 3 + 1 + 2 + 9 = 19 components. Store them in one `samurai` vector field. Define the
component offsets in one place (`layout.hpp`); never hard-code an index elsewhere.

In the papers' notation, a^b = A[b][0], b^b = A[b][1], c^b = A[b][2].

### 4.2 Constitutive relations

Stiffened gas for each phase:

    rho_l e^h_l = (p_l + gamma_l p_inf_l) / (gamma_l - 1)
    c_l^2 = gamma_l (p_l + p_inf_l) / rho_l

Shear part for a solid phase (Ndanou 2015 eq. 5):

    G = A^T A,   g = G / |det A|^(2/3),   j2 = tr(g^2)
    e^e = mu / (8 rho0) (j2 - 3)
    S   = -(mu rho_l / (2 rho0)) (g^2 - (j2/3) I)        (deviatoric, symmetric)

Mixture:

    rho = sum (alpha rho)_l,   Y_l = (alpha rho)_l / rho
    sigma = -(sum alpha_l p_l) I + chi * sum_s alpha_s S_s
    rho E = sum (alpha rho e^h)_l + chi * sum_s (alpha rho)_s e^e_s + rho |v|^2 / 2

chi is the fluidisation flag of section 4.4 (1 in solid cells, 0 in fluid cells).

Favrie 2009 uses a slightly different shear energy, mu/(4 rho0) tr((g - I)^2). Both give Hooke's law at small strain
and differ by about 1 % in shear stress at 15 % strain. Keep the shear energy in one function so the other can be
added later.

### 4.3 Wave speeds

For a face of normal direction d, the acoustic tensor (per unit mass) is

    K = c_f^2 e_d e_d^T + chi * sum_s Y_s D_s,        c_f^2 = sum Y_l c_l^2
    D_s = mu/(2 rho0) [ g_dd g + (g e_d)(g e_d)^T + g^2
                        - ((g^2 e_d) e_d^T + e_d (g^2 e_d)^T)/3 + (j2/9) e_d e_d^T ]

At g = I this gives 4 mu/(3 rho0) along d and mu/rho0 across, as expected. It was derived for this plan (not in the
papers) and checked against a direct uniaxial calculation.

    c_long = sqrt( max_i sum_j |K_ij| )                     (Gershgorin bound, safe for HLLC)
    c_t    = chi * sqrt( max(K_t1t1, K_t2t2) + |K_t1t2| )   (t1, t2: the two directions other than d)

Do not use the small-strain formula alone: at 15 % compression it underestimates the elastic part by a factor 1.6.

### 4.4 Slip treatment: fluidisation

Problem: a residual solid fraction gives the gas shear stiffness, and a 3-wave HLLC couples transverse velocity
across contacts. Favrie 2009 and 2012 fix this with a level set and a ghost-fluid step. We use a local rule instead:

    chi = 1 if (sum of solid volume fractions) >= alpha_lim, else 0

In a cell with chi = 0 the solid has no shear strength: its stress and shear energy are not counted, its shear
impedance is zero, and at every relaxation step each A_s is reset to `|det A_s|^(1/3) I`. The elastic energy removed
becomes heat through the energy reset.

Combined with the shear solver of section 5.3, a face between a solid cell and a fluid cell gets zero shear traction
and the free-surface tangential velocity, with one flux per face (conservative).

Default `alpha_lim = 1e-2`, as a run-time parameter. Prototype evidence (1000 cells, order 2):

| Test | alpha_lim | Result |
|---|---|---|
| Favrie 6.1, residual 1e-7 | no treatment | spurious shear strain 1.4, v overshoots to 245 m/s (should be 100), p from 0.75e5 to 1.08e5 |
| same | 0.5 | strain 2.7e-7, shear stress 2.5e4 Pa, v exact |
| same | 1e-2 | strain 1.7e-6, shear stress 7.5e4 Pa |
| same | 1e-4 | strain 6e-4 |
| Favrie 6.2 at 180 microseconds | 0.5 | 27 kg/m2 of solid ejected from the free surfaces |
| same | 1e-2 | 0.39 kg/m2 |
| same | 1e-3 | 0.045 kg/m2 |

A large threshold dumps the elastic energy of strongly sheared interface cells into heat and blows them off; a small
one lets gas-dominated cells keep a strain. 1e-2 to 1e-3 is the compromise. Fallback if 2D shows problems: the
ghost-fluid step of Favrie 2012 section 5.2.

## 5. Numerical method

### 5.1 Overview of one time step (SSP-RK2)

    dt = cfl * h_min / max over cells and directions of (|v_d| + c_long)        cfl = 0.4
    U1 = relax( U + dt L(U) )
    U  = relax( (U + U1 + dt L(U1)) / 2 )

L is the sum over mesh directions of the flux operator of sections 5.2 to 5.4. `relax` is section 5.5. First order
(`--order 1`) is forward Euler with zero slopes; keep it as an option for debugging.

The papers use MUSCL-Hancock. SSP-RK2 with relaxation after each stage is used here because it maps directly onto
Samurai's scheme operators. This choice is validated in the prototype.

### 5.2 Reconstruction

Stencil of four cells per face (j-1, j, j+1, j+2), ghost width 2. Convert the four states to primitive variables
`(alpha_l, rho_l, v, p_l, A_s)`, compute minmod slopes in cells j and j+1, and build the left state `W_j + s_j/2` and
the right state `W_j+1 - s_j+1/2`. Then divide the reconstructed alpha_l by their sum. chi is not reconstructed: the
left state uses chi of cell j, the right state chi of cell j+1.

Favrie 2009 (appendix B) sets the slopes of the cobasis to zero. The prototype reconstructs it with minmod without
trouble. Provide `--muscl-cobasis 0|1`, default 1; set 0 if interface problems appear in 2D.

### 5.3 Face solver for direction d

Notation: n = v_d, K in {L, R}, sigma_n = sigma_dd, tau_t = sigma_td for the two tangential directions t.

Longitudinal part (HLLC, tangential momentum only advected):

    S_L = min(n_L - c_long_L, n_R - c_long_R),   S_R = max(n_L + c_long_L, n_R + c_long_R)
    S_M = [sigma_n_L - sigma_n_R + rho_L n_L (S_L - n_L) - rho_R n_R (S_R - n_R)]
          / [rho_L (S_L - n_L) - rho_R (S_R - n_R)]
    upwind side K = L if S_M >= 0 else R;  supersonic if (K = L and S_L >= 0) or (K = R and S_R <= 0)
    fac = (S_K - n_K)/(S_K - S_M)   (1 if supersonic),   n* = S_M   (n_K if supersonic)

    conservative fluxes, F = F_K + S_K (U*_K - U_K)   (F = F_K if supersonic):
      (alpha rho)_l :  F_K = (alpha rho)_l n,            U* = (alpha rho)_l fac
      rho v_d       :  F_K = rho n^2 - sigma_n,          U* = rho fac S_M
      rho v_t       :  F_K = rho n v_t,                  U* = rho fac v_t
      rho E         :  F_K = (rho E - sigma_n) n,
                       U* = fac [rho E + (S_M - n)(rho S_M - sigma_n/(S_K - n))]

    advected quantities, flux = n* phi*_K:
      alpha_l          : phi* = alpha_l
      A, column d      : phi* = A[.][d] fac          (behaves like a density)
      A, columns t     : phi* = A[.][t]
      (alpha rho e^h)_l: phi* = alpha_l (p*_l + gamma_l p_inf_l)/(gamma_l - 1),  rho*_l = rho_l fac,
                         p*_l + p_inf = (p_l + p_inf) [(g-1) rho_l - (g+1) rho*_l] / [(g-1) rho*_l - (g+1) rho_l]
                         (phase Hugoniot; if the denominator is not negative or the result not positive,
                          use the isentrope (p_l + p_inf) fac^gamma)

Shear part, for each tangential direction t, with Z_K = rho_K c_t_K:

    if Z_L + Z_R > 0:
        v_t*   = (Z_L v_t_L + Z_R v_t_R + tau_t_R - tau_t_L) / (Z_L + Z_R)
        tau_t* = (Z_R tau_t_L + Z_L tau_t_R + Z_L Z_R (v_t_R - v_t_L)) / (Z_L + Z_R)
    else:
        v_t* = (v_t_L + v_t_R)/2,   tau_t* = 0
    add to the fluxes:   rho v_t : -tau_t*        rho E : -tau_t* v_t*

The face returns the fluxes F and the face velocity vector `v* = (n* along d, v_t* along the tangents)`.

Summing the two parts unsplit is a design choice of this plan; Favrie 2014 applies the sub-systems one after the
other. The unsplit sum is validated in the prototype at CFL 0.4.

### 5.4 Cell update and Samurai mapping

For a cell i and direction d, with faces i-1/2 and i+1/2:

    dU_i/dt = -(F_i+1/2 - F_i-1/2)/h - C_i (v*_i+1/2 - v*_i-1/2)/h

where the non-conservative coefficient C_i uses the values of cell i:

| Quantity | Extra term in dU/dt |
|---|---|
| alpha_l | `+ alpha_l * dn` |
| (alpha rho e^h)_l | `- alpha_l p_l * dn` |
| A, column t (t != d) | `+ A[.][t] * dn` |
| A, column d | `- sum over t of A[.][t] * dv_t` |

with `dn = (n*_i+1/2 - n*_i-1/2)/h` and `dv_t = (v_t*_i+1/2 - v_t*_i-1/2)/h`.

Samurai: use a non-linear flux-based scheme with the non-conservative `flux_function`, which fills a
`samurai::FluxValuePair<cfg>`: element 0 is the contribution of the face to its left cell, element 1 to its right
cell, and the scheme output is a divergence (the caller writes `U - dt * scheme(U)`). For a face between cells L
and R of the stencil:

    flux[0] =  (F + C_L applied to v*)       C_L from stencil cell j
    flux[1] = -(F + C_R applied to v*)       C_R from stencil cell j+1

where "C applied to v*" means the terms of the table with the differences replaced by the face value, with the
opposite sign of the table (the table gives dU/dt, the scheme gives the divergence). For conservative components
this reduces to flux[1] = -flux[0].

This sign convention was read in `flux_definition.hpp` and `flux_based_scheme__nonlin.hpp` (samurai 0.30.1) but not
run. First task of milestone M2: confirm it with the uniform-advection test, which fails visibly if a sign is wrong.

At a level jump Samurai evaluates the flux function on fine cells, some of which are ghosts filled by prediction;
C then comes from a ghost. This is consistent to second order and needs no special code.

### 5.5 Relaxation step (cell by cell)

1. Normalise the alpha_l to sum 1. Compute chi. For chi = 0, reset each `A_s` to `|det A_s|^(1/3) I`.
   (Plasticity will be inserted here later.)
2. Pressure relaxation. With p0_l the phase pressures (clipped to stay above -p_inf_l), solve for p

       f(p) = sum_l alpha_l (p0_l + gamma_l p_inf_l + (gamma_l - 1) p) / (gamma_l (p + p_inf_l)) - 1 = 0

   f is decreasing and convex on p > p_min = max_l(-p_inf_l). Newton from `max(sum alpha_l p0_l, p_min + small)`;
   if an iterate leaves the domain, take the midpoint with p_min. Then multiply each alpha_l by its term of the sum.
   Starting far left of the root Newton only doubles p per iteration (up to about 40 iterations when the gas is
   near vacuum); for two phases the equation is a quadratic and may be solved in closed form.
3. Energy reset. `rho e^h = rho E - rho |v|^2/2 - chi * sum (alpha rho)_s e^e_s`, then

       p_e = (rho e^h - sum alpha_l gamma_l p_inf_l/(gamma_l - 1)) / (sum alpha_l/(gamma_l - 1))

   If p_e > p_min, set every `(alpha rho e^h)_l = alpha_l (p_e + gamma_l p_inf_l)/(gamma_l - 1)`. Otherwise keep the
   relaxed pressure p of step 2 for the phase energies and recompute rho E from them (local loss of energy
   conservation). Count these cells and print the count in the log.

Step 3's fallback matters when the gas is nearly a vacuum (void opening with a residual of 1e-7): energy errors that
are negligible for the solid make p_e negative. Without the fallback both phases get a negative pressure and the
void stops growing. With a residual of 1e-4 the prototype never used it in the cavitation test.

### 5.6 Boundary conditions

Transmissive: copy the boundary cell into both ghost layers (copy `TransmissiveImpl` of `bc_euler.hpp` with R = 2;
it is generic in the number of components). Symmetry plane or rigid wall of normal d, for later (half-domain runs,
Taylor impact): mirror the cells, flip the sign of v_d, of A[d][i] for i != d and of A[b][d] for b != d.

## 6. Code organisation

    CMakeLists.txt            find_package(samurai), C++20, Release by default, one executable `hvi`
    README.md
    src/layout.hpp            component offsets, N, Ns, material table
    src/physics.hpp           header-only kernel, no Samurai include: Material, 3x3 helpers, prim <-> cons,
                              elastic(), acoustic tensor, wave speeds, face solver (direction d), relax()
    src/scheme.hpp            Samurai flux-based scheme: stencil -> primitives -> MUSCL -> face solver -> FluxValuePair
    src/bc.hpp                transmissive (and later symmetry) conditions, R = 2
    src/cases.hpp             initial conditions of all tests
    src/main.cpp              options, mesh, time loop, adaptation, outputs
    tests/test_physics.cpp    kernel checks against numbers printed by the prototype
    scripts/plot_1d.py, scripts/schlieren.py, scripts/compare_proto.py
    prototype/hvi1d.py        reference (exists)

`physics.hpp` must be templated on nothing mesh-related and use plain `std::array<double, ...>`, so that it can be
tested without Samurai and reused in 1D, 2D and 3D. The mesh dimension `dim` is a compile-time constant of
`main.cpp` (CMake option `HVI_DIM`, executables `hvi1d`, `hvi2d`).

Build, as in the Euler example (conda env `cwisam_catapult`, samurai 0.30.1):

    conda activate cwisam_catapult
    cmake -S . -B build -DCMAKE_PREFIX_PATH=$CONDA_PREFIX && make -C build -j

Samurai patterns (all present in the Euler example):

    auto config = samurai::mesh_config<dim>().min_level(lmin).max_level(lmax).max_stencil_radius(2);
    config.scaling_factor(base);                     // cell size = base / 2^level
    auto mesh = samurai::mra::make_mesh(samurai::Box<double, dim>(lo, hi), config);
    auto U = samurai::make_vector_field<double, ncomp>("U", mesh);
    samurai::make_bc<hvi::Transmissive>(U);   U1.copy_bc_from(U);
    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, 4, Field, Field>;
    samurai::FluxDefinition<cfg> def;  def[d].flux_function = [...](samurai::FluxValuePair<cfg>& flux,
        const samurai::StencilData<cfg>& data, const samurai::StencilValues<cfg>& u) {...};
    auto L = samurai::make_flux_based_scheme(def);   // ghosts are updated automatically when it is applied

Mesh sizes: with `base = 0.05` the interfaces of every 1D test and of the 2D case fall on cell faces at every level
(20 x 2^level cells per metre: level 6 is 1280, level 9 is 10240). For 2D use `base = 0.05` as well (level 5 is
640^2, level 6 is 1280^2).

Adaptation (pattern of `turek_fsi_fluid.cpp`): adapt on a normalised indicator field and carry U along.

    ind = ( sum of solid alpha,  log(rho)/log(rho_s0/rho_g0),  v / U_ref,  p / (rho_s0 c_s0 U_ref) )
    auto mra_cfg = samurai::mra_config().epsilon(1e-3).regularity(1.);
    auto mr_adapt = samurai::make_MRAdapt(ind);
    mr_adapt(mra_cfg, U);  then resize the scratch fields, then apply relax() to all cells

Adapt once per time step, before the step. The density and volume-fraction jumps keep the interfaces at the finest
level. The call to relax() after adaptation restores pressure equilibrium in cells created by prediction.

Outputs: as in `turek_fsi_fluid.cpp`, a `.npy` list of cells `(x, [y,] h, alpha_s, rho, u, v, w, p, sigma components,
A components, level)` per output time, plus a log line per step (t, dt, cells, total mass, momentum, energy, number
of floored cells). `scripts/schlieren.py` rasterises the cells on the finest grid and plots `exp(-k |grad rho| /
max|grad rho|)`.

Machine: 8 cores, 7 GB of RAM. A uniform 1000^2 run of case A is of the order of 10^4 steps and many hours on one
core; develop at levels 4 to 5 (320^2, 640^2) and keep the finest runs for the adapted mesh.

## 7. Verification ladder (1D, copper and air)

All are implemented in the prototype (`python prototype/hvi1d.py <case> --plot out.png`). The first four have been
run; the reference numbers are for 1000 cells, order 2, CFL 0.4, alpha_lim as stated.

| # | Case (prototype name) | Setup | Acceptance |
|---|---|---|---|
| 1 | `advect` (Favrie 5.1) | interface at 0.5, u = 1000, p = 1e5, residual 1e-7, t = 150 microseconds | u and p uniform to round-off, at order 1 and 2 |
| 2 | `slip-advect` (Favrie 6.1) | as 1, v = 100 in the solid, 0 in the air | v within [0, 100.1]; max shear stress below 1e5 Pa; p within [1e5, 1.02e5] (alpha_lim 1e-2) |
| 3 | `shear-slab` (Favrie 6.2) | slab [0.15, 0.85], v = 1000 on [0.35, 0.65], residual 1e-7 | at 54.6 microseconds: plateau v = 507, A[1][0] = -/+0.157, shear stress 14.65 GPa, fronts at 0.175 and 0.525 (and mirrored); L1 error on v against the linear solution about 10 m/s |
| 4 | `double-expansion` (Favrie 5.4) | residual 1e-4, u = -3000 left, +3000 right, 2000 cells, t = 88 microseconds | a void opens at the centre (alpha_s min about 2e-3, rho min about 13), pressure positive everywhere |
| 5 | `solid-shear-tube` (Favrie 5.5) | pure copper, p = 1e8 / 1e5, v = 0 / 100 | five waves, matches the prototype |
| 6 | `solid-gas-tube` (Favrie 5.3) | copper at 5e9 Pa / air at 50 kg/m3, interface 0.6 | runs to 87 microseconds, matches the prototype |
| 7 | `impact1d` | copper slab at 800 m/s onto a copper slab at rest, in air, residual 1e-4 | matches the prototype; contact velocity 400 m/s |

Cases 5 to 7 exist in the prototype but have not been run yet: run them in the prototype first and fix the prototype
if needed before using it as a reference.

About case 3: at 15 % shear this model is strongly nonlinear. The sheared zones want a negative pressure of about
-1.5 GPa (normal-stress effect of finite shear), the solid cavitates to alpha_s of about 0.96 by 180 microseconds
and longitudinal waves of +/-90 m/s appear. The result therefore does not look like the linear curves of Fig. 14 at
the later time. Use the early time for validation, and optionally v = 10 m/s for a clean comparison with the linear
solution (`dalembert_slab` in the prototype).

Each C++ case is compared with the prototype by `scripts/compare_proto.py` (same cells, same order): relative L1
difference below 1e-3 on rho, u, v, p. Differences at round-off level are expected only at order 1.

## 8. Milestones

Each milestone ends with its tests passing and a short entry in `README.md`.

- **M0, scaffold.** CMake project, empty executable linking Samurai, builds in `cwisam_catapult`.
- **M1, kernel.** `physics.hpp` and `tests/test_physics.cpp`. Checks: prim/cons round trip; S = 0 and e^e = 0 at
  A = I; simple shear `A[1][0] = -s` gives `S_01 = mu s (1 + s^2/2)`; acoustic tensor at A = I; uniaxial value
  `D_dd = (mu/rho0)(11/9 a^(8/3) + 1/9 a^(-4/3))` for A = diag(a, 1, 1); the face solver returns the exact flux for
  equal states; relax() leaves an equilibrium state unchanged; a few face-solver outputs equal to the prototype's.
- **M2, 1D uniform.** `scheme.hpp`, `main.cpp` with `min_level = max_level`. First confirm the FluxValuePair sign
  convention with case 1. Then cases 1 to 7 against the prototype, at order 1 and 2.
- **M3, 1D adapted.** Multiresolution on cases 2, 3, 4, 7. Acceptance: same results as the uniform run at the finest
  level within a few percent in L1, interfaces at the finest level, no negative density or volume fraction.
- **M4, 2D uniform.** Each 1D case run along x and along y gives the same profile, with the transverse components
  exactly zero where they should be. Then the impact, case A, at level 4 (320^2) to 100 microseconds: check total
  mass and x-momentum conservation (momentum can leave only through the boundaries) and symmetry about y = 0.5.
- **M5, target.** Impact with adaptation, cases A, B, C, to the last paper time, schlieren images at the paper
  times, comparison with a uniform run at the same finest level, cell counts and run times in the README.
- **M6, plasticity.** Maxwell relaxation to the von Mises surface (Favrie 2012 sections 3.4 and 5.3; singular value
  decomposition of each A_s), inserted in step 1 of relax(). Tests: Wilkins flying plate (Favrie 2012 section 6.1),
  then the impact of Favrie 2012 section 6.2.
- **M7, later.** Second solid phase (Ndanou 2015 section 6.3), symmetry boundary, 3D, MPI.

## 9. Risks and how to handle them

| Risk | Sign | Action |
|---|---|---|
| Wrong sign in the non-conservative pair | case 1 loses uniform u, p | fix before anything else |
| Multiresolution prediction overshoots near interfaces | alpha outside [0, 1], negative phase density after adaptation | interfaces must be at the finest level (check the indicator); clip alpha to [1e-12, 1] and renormalise before relax(); as a last resort lower epsilon |
| Negative energy-based pressure | floored-cell count grows | expected in void opening; if it happens in compressed solid, reduce CFL and check the wave-speed bound |
| Elastic energy dumped at fluidised cells | spray of solid leaving the surfaces, hot gas pockets | lower alpha_lim (1e-3); compare with alpha_lim = 0.5 |
| Shear waves in the air | transverse velocity spreading in the gas at km/s | chi not applied somewhere (stress, energy, impedance, reset) |
| Case C (mu = 0) | division by zero in shear solver | handled by the Z_L + Z_R = 0 branch; test it in M1 |
| Corners and oblique interfaces | staircase artefacts, filaments (also seen in the paper, Fig. 17) | accept; document |
| Cost | hours per run | develop at 320^2; profile the stencil-to-primitive conversion first |

Decisions taken in this plan that differ from the papers, all switchable: unsplit longitudinal + shear flux instead
of the 3-wave HLLC; fluidisation instead of level set and ghost fluid; SSP-RK2 instead of MUSCL-Hancock; limited
slopes on the cobasis; relaxed-pressure fallback in the energy reset; Ndanou's shear energy for a 2009 test.
