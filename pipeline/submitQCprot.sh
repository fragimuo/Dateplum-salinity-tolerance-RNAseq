#!/bin/bash -l

set -ex

proj=u2018012
mail="francisco.gil.munoz@slu.se"
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/functional
trans=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity
index=/mnt/picea/storage/reference/UniRef90/201701/indices/blast+2.6.0/uniref90

mkdir -p $out

sbatch -A $proj --mail-user $mail -e $out/blastprot.err -o $out/blastprot.out \
$UPSCb/pipeline/runQCprot.sh $trans \
$index $out