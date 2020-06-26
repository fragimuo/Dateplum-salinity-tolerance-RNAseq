#!/bin/bash -l

set -ex

proj=u2018012
mail=francisco.gil.munoz@slu.se

in=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trimmomatic
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/fastqc/trimmomatic

module load bioinfo-tools FastQC

mkdir -p $out

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*.fq.gz"`; 
do
  fnam=$(basename ${f/.fq.gz/})
  sbatch -A $proj --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out $UPSCb/pipeline/runFastQC.sh $out $f
done
