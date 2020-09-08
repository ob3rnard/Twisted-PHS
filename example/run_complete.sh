#!/bin/bash

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <nf> <bsz> <N> , where :\n\t<nf> is field label (eg. 'z23' or 'n37'),\n\t<bsz> is the challenge prime norm bit size,\n\tN is the number of challenges."
    exit 2
fi

nf=$1;
bsz=$2;
N=$3;

expl_folder=`dirname $0`;

# setting
cd ${expl_folder}/../scripts/ ;
mkdir -p ../data/$nf/ ;

# places
./nf_get_places.sage $nf ;
./nf_get_chals.sage  $nf $bsz $N ;

# S-units, ClDL for phs/opt/tw
magma -b nf:=$nf sunits.magma ;
for met in tw opt phs; do
    magma -b nf:=$nf typ:=$met chal:="../data/$nf/$nf.chal_b${bsz}_n${N}" cldl.magma ;
done;

# Compute Log-sunit lattices, reduce
./lat_logsunit.sage $nf ;
./lat_reduce.sage $nf ;

# Geometry can be evaluated here
./eval_geo.sage $nf ;
./eval_gs_norms.sage $nf ;
./plot_gs_norms.sh $nf ;

# Solving idSVP (N targets*~2s* 5 trials*5 methods --> ~10 min
echo "This will take up to 10 minutes, but please look at the previous graphs at ../data/$nf/*.png";
./solve_idsvp.sage $nf $bsz $N ;



exit 0;
