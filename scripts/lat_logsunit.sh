#!/bin/bash


for nf in "$@"; do
    echo "'$nf' field: compute log s-unit lat 'tw/opt/phs/opt0/phs0'";
    ./lat_logsunit.sage $nf 2>&1 1> ../logs/$nf.latlog &
done;
    
exit 0;
