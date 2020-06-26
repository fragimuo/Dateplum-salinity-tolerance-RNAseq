#!/bin/bash -l

## error and verbose
set -ex

module load R

proj=u2018012
mail="francisco.gil.munoz@slu.se"
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/hmmer
tmp=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/hmmer/tmp/

if [ ! -d $out ]; then
    mkdir -p $out
fi

if [ ! -d $tmp ]; then
    mkdir -p $tmp
fi

sbatch -A $proj --mail-user $mail -e $out/PrepprotB2GO.err -o $out/PrepprotB2GO.out \
-J PrepprotB2GO $UPSCb/projects/persimmon/pipeline/prepareprotBLAST2GO.R