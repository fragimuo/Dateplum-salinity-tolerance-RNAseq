#!/bin/bash -l

sbatch -e /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/HMMER.err \
-o /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/HMMER.out \
$UPSCb/projects/persimmon/pipeline/runHMMER.sh /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/ /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity/Trinity.format.fasta.transdecoder.pep