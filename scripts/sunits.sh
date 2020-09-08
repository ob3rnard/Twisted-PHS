#!/bin/bash

for nf in "$@"; do
    echo "Launching S-Units for $nf";
    magma -b nf:=$nf sunits.magma > ../logs/${nf}.sulog &
done;

exit 0;
