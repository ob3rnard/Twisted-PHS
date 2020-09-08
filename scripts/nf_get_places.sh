#!/bin/bash


for nf in "$@"; do
    echo "'$nf' field: get places for 'tw/opt/phs'";
    ./nf_get_places.sage $nf  2>&1 1> ../logs/$nf.fblog &
done;

exit 0;
