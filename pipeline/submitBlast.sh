#!/bin/bash -l

set -ex

proj=u2018012
mail="francisco.gil.munoz@slu.se"
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/functional
trans=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity/Trinity.format.fasta
index=/mnt/picea/storage/reference/NCBI/20171103/nr

mkdir -p $out

sbatch -A $proj --mail-user $mail -e $out/blast.err -o $out/blast.out \
$UPSCb/projects/persimmon/pipeline/runBlastPlus.sh -p 8 blastx $trans \
$index $out