#!/bin/bash -l

sbatch -e /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/trinotate.err \
-o /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/trinotate.out \
$UPSCb/projects/persimmon/pipeline/runTrinotate.sh /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity/Trinity.format.fasta