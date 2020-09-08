#!/bin/bash

# At least one lattice in (only tw/opt0/phs0 for z61)
# n23 n29 n31 n37 n41 n43 n47   z23 z29 z31 z37 z41 z43 z47 z53 z59 z61
for nf in "$@"; do
    echo "'$nf' field: reduce log s-unit lat 'tw/opt/phs/opt0/phs0'";
    ./lat_reduce.sage $nf 2>&1 1> ../logs/$nf.bkzlog &
done;

exit 0;
