#!/bin/bash

bsz=$1;
N=$2;

for nf in "${@:3}"; do
    echo "'$nf' field: get places for 'tw/opt/phs'";
    ./nf_get_chals.sage $nf $bsz $N  2>&1 1> ../logs/$nf.challog_${bsz}_${N} &
done;

exit 0;
