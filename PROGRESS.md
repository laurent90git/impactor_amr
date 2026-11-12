# Progress log (night of 2026-11-12)

Resume point for any restart: read this file, then DEVELOPMENT_PLAN.md.

## Done
- M0/M1: C++ code in `src/` (physics kernel, numerics, samurai glue, cases, driver); `tests/test_physics.cpp` passes (12 checks).
- M2: 1D uniform. All 7 prototype cases match `prototype/hvi1d.py` to round-off (rel. L1 diff 1e-9 .. 1e-14) at 320 cells
  (`scripts/compare_proto.py <case> out/t/<case>_000.npy --n 320`). This also confirms the samurai flux-pair signs.

## Rules for running on this laptop
- 8 cores, 7 GB RAM: `ulimit -v 3000000; nice ...`, one core per run during development. Parallel runs allowed later for the final cases.
- M3: 1D multiresolution: 9-21x fewer cells, 0.1-3 % L1 from uniform finest level.
- Plasticity check (Ndanou 2015 6.2.1, case `kanel1d`): precursors within 0.2 %, plastic plateau within 0.4 % of the paper's analytical values.
- M4: 2D builds; 1D cases along x and y reproduce 1D; 2D impact: mass / x-momentum exact, energy +1e-5 (cavitation floor).
- Docs: Sphinx sources in docs/ (build: `~/miniconda3/envs/catapult_doc/bin/sphinx-build -b html docs docs/_build/html`).

## Running (started ~01:00)
- `scripts/run_production.sh` (5 runs, out/prod/*.txt): impact2009 A/B/C, impact2012 (plastic), penetration2015. Level 5, MR from level 2.
- `scripts/make_doc_figures_1d.py` -> docs/_static/fig (1D figures and tables).

## Next
- 2D figures for docs/impacts.rst from out/prod (scripts/schlieren.py), compare with the papers' figures.
- If a run fails: look at out/prod/<X>.txt, rerun that line of run_production.sh.
