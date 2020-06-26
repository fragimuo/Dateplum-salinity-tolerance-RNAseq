#!/bin/bash -l

## error and verbose
set -ex

module load R

proj=u2018012
mail="francisco.gil.munoz@slu.se"
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/p/
tmp=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/ptmp/

if [ ! -d $out ]; then
    mkdir -p $out
fi

if [ ! -d $tmp ]; then
    mkdir -p $tmp
fi

sbatch -A $proj --mail-user $mail -e $out/PrepB2GO.err -o $out/PrepB2GO.out \
-J PrepB2GO $UPSCb/projects/persimmon/pipeline/prepareBLAST2GO.R