#!/bin/bash
# Large 3D run (target ~1M cells): elastic-plastic disc plate, base 0.04 m, levels 1..3 (finest cell 5 mm),
# OpenMP build (build-omp, -DSAMURAI_WITH_OPENMP=ON), 6 threads, memory capped at 6 GB, checkpoint every 15 min.
cd "$(dirname "$0")/.."
ulimit -v 6000000
mkdir -p out/3dbig
OMP_NUM_THREADS=6 nice -n 10 ./bin/hvi3d_omp --case impact3d-disc-plastic --base 0.04 --level 3 --coarse-level 1 \
  --tf 300e-6 --nout 30 --log-every 50 --checkpoint-every 900 --out out/3dbig "$@" > out/3dbig/run.txt 2>&1
