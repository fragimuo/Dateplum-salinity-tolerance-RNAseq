#!/bin/bash -l

set -ex

proj=u2018012
mail="francisco.gil.munoz@slu.se"
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity
in=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity
trans=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity/Trinity.fasta

sbatch -A $proj --mail-user $mail -e $out/ExN50.err -o $out/ExN50.out \
$UPSCb/projects/persimmon/pipeline/runExN50.sh $in/Trinity_trans.isoform.TMM.EXPR.matrix \
$trans $out