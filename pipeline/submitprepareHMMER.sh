#!/bin/bash -l

sbatch -e /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/prepareHMMER.err \
-o /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/prepareHMMER.out \
$UPSCb/projects/persimmon/pipeline/runprepareHMMER.sh 