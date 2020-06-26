#!/bin/bash -l

set -ex

proj=u2018012
mail="francisco.gil.munoz@slu.se"
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/Diff
samples=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity/samples.txt
trinity_index=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity
  

if [ ! -d $out ]; then
    mkdir -p $out
fi
  
sbatch -A $proj --mail-user $mail -e $out/DESeq.err -o $out/DESeq.out \
-J DESeq $UPSCb/projects/persimmon/pipeline/runDESeq.sh $trinity_index $samples $out

