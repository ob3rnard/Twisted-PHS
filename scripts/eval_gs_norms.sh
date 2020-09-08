#!/bin/bash

# At least one lattice in (only tw/opt0/phs0 for z61)
# n23 n29 n31 n37 n41 n43 n47   z23 z29 z31 z37 z41 z43 z47 z53 z59 z61
for nf in "$@"; do
    echo "'$nf' field: compute gs-log norms of lats 'tw/opt/phs/opt0/phs0'";
    ./eval_gs_norms.sage $nf 2>&1 1> ../logs/$nf.gsnlog &
done;

exit 0;
