#!/bin/bash
# Production runs (2D, multiresolution levels 2..5), one core each, niced, memory-capped.
cd "$(dirname "$0")/.."
ulimit -v 2500000
B=./build/hvi2d
O=out/prod
nice -n 10 $B --case impact2009 --mu 9.2e10 --tag A --level 5 --coarse-level 2 --paper-times --log-every 500 --out $O > $O/A.txt 2>&1 &
nice -n 10 $B --case impact2009 --mu 1e9 --tag B --level 5 --coarse-level 2 --paper-times --log-every 500 --out $O > $O/B.txt 2>&1 &
nice -n 10 $B --case impact2009 --mu 0 --tag C --level 5 --coarse-level 2 --paper-times --log-every 500 --out $O > $O/C.txt 2>&1 &
nice -n 10 $B --case impact2012 --tag plastic --level 5 --coarse-level 2 --paper-times --log-every 500 --out $O > $O/P.txt 2>&1 &
nice -n 10 $B --case penetration2015 --level 5 --coarse-level 2 --tf 120e-6 --nout 6 --log-every 500 --out $O > $O/N.txt 2>&1 &
wait
