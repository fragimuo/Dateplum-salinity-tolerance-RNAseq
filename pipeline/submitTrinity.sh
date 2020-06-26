#!/bin/bash -l

sbatch -e /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity.err \
-o /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity.out \
$UPSCb/projects/persimmon/pipeline/runTrinity.sh