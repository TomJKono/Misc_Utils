#!/bin/bash

#   Thomas Kono
#   Saint Paul, MN
#   2015-09-09

#   A script to delete all owned jobs from the University of Minnesota
#   Supercomputing Institute (MSI) queue. Designed to work with the portable
#   batch system (PBS).

set -e
set -u
set -o pipefail

USERNAME=$(whoami)

ALLJOBS=($(qstat -a -u ${USERNAME} | tail -n +6 | cut -f 1 -d ' ' | cut -f 1 -d '.'))

echo "Will delete the following jobs:"
echo ${ALLJOBS[@]}

for j in ${ALLJOBS[@]}
    do qdel $j
done
