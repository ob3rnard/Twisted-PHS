#!/bin/bash

bsz=$1;
N=$2;

for nf in "${@:3}"; do
    echo "Launch ClDL for $nf (bsz=$bsz,#=$n) [3 thr.]";
    for typ in tw opt phs; do
        magma -b nf:=$nf typ:=$typ chal:="../data/$nf/$nf.chal_b${bsz}_n$N" cldl.magma > ../logs/${nf}_$typ.cldllog_${bsz}_${N} &
    done;
done;

exit 0;
