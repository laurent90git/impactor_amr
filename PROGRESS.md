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

## Night, ~02:00: time-step collapse in B, impact2012, penetration2015 (fixed in source, being re-run)
- Diagnosis (case B at 67.8 us): copper cell at alpha = 1 drained to rho = 34 kg/m3 at p = 0, air crushed; the energy
  floor injected energy -> hot light copper with c = 65 km/s -> dt -> 0.
- Fix: relaxation step 5 opens a void (compacts the stiff phases, gas fills) instead of injecting energy (`voided`
  counter). Also: slip switch per solid phase (needed with 2 solids). 1D results unchanged.
- Runs A and C (old binary in build/) continue; B and 2012 re-run with build-dev/hvi2d into out/dbg.

## Morning 2026-11-13
- Restart capability: `--checkpoint-every <s>` (default 900), checkpoint + clean stop on SIGTERM/Ctrl-C, `--resume <prefix>_ckpt`.
  Bitwise identical to uninterrupted runs (1D and 2D, adaptive). Binaries bin/hvi*_v5.
- PAUSED (SIGSTOP) at the user's request: pid 149004 (impact2012 P4 level 5, ~612/1100 us) and 151651 + wrapper 151650
  (penetration2015 GIF run, ~80/100 us). These are v4 binaries without checkpoints: resume with
  `kill -CONT 149004 151650 151651`.
