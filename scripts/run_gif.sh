#!/bin/bash
# Level-4 runs with 40 evenly spaced outputs, for the animations (scripts/make_gif.py). Two at a time.
cd "$(dirname "$0")/.."
ulimit -v 2500000
B=./bin/hvi2d_v4
O=out/gif
cat <<LIST | xargs -P 2 -I{} bash -c "nice -n 10 $B {} --level 4 --coarse-level 2 --nout 40 --log-every 1000 --out $O > $O/\$(echo {} | awk '{print \$NF}').txt 2>&1"
--case impact2009 --mu 9.2e10 --tf 609e-6 --tag A
--case impact2009 --mu 0 --tf 1.03e-3 --tag C
--case impact2009 --mu 1e9 --tf 710e-6 --tag B
--case impact2012 --tf 1.1e-3 --tag P
--case penetration2015 --tf 100e-6 --tag N
LIST
