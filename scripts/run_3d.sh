#!/bin/bash
# Coarse 3D impacts (multiresolution levels 1..2, finest cell 12.5 mm), one after the other on one core.
cd "$(dirname "$0")/.."
ulimit -v 2500000
for c in impact3d-disc impact3d-square impact3d-disc-plastic; do
  nice -n 10 ./bin/hvi3d_v6 --case $c --level 2 --coarse-level 1 --tf 300e-6 --nout 30 --log-every 200 --out out/3d > out/3d/$c.txt 2>&1
done
