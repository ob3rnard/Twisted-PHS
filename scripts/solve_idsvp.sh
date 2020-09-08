#!/bin/bash

bsz=$1
N=$2

# cldl sol for 
# n23 n29 n31 n37   z23 z29 z31 z37 z41 z43 z47 #z53
for nf in "${@:3}"; do
    echo "Field:'$nf', solving idsvp for 'tw/opt/phs/opt0/phs0'";
    ./solve_idsvp.sage $nf $bsz $N 2>&1 1> ../logs/$nf.sol_b${bsz}_n${N} &
done;

exit 0;
