# Progress log (night of 2026-11-12)

Resume point for any restart: read this file, then DEVELOPMENT_PLAN.md.

## Done
- M0/M1: C++ code in `src/` (physics kernel, numerics, samurai glue, cases, driver); `tests/test_physics.cpp` passes (12 checks).
- M2: 1D uniform. All 7 prototype cases match `prototype/hvi1d.py` to round-off (rel. L1 diff 1e-9 .. 1e-14) at 320 cells
  (`scripts/compare_proto.py <case> out/t/<case>_000.npy --n 320`). This also confirms the samurai flux-pair signs.

## Rules for running on this laptop
- 8 cores, 7 GB RAM: `ulimit -v 3000000; nice ...`, one core per run during development. Parallel runs allowed later for the final cases.
