#!/bin/bash -l

## stop on error
set -e
# use set -ex for debugging instead

## args
mail="francisco.gil.munoz@slu.se"
proj=u2018012
  

## create the in and out dir
in=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/raw
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trimmomatic
mkdir -p $out

module load bioinfo-tools
module load Trimmomatic 

## select all files
i=0
id=422141
for f in `find $in -name "*.fq.gz"`; do echo "$(basename ${f//_[1,2].fq.gz/})" ; done | sort | uniq | while read line;
do
i=$(expr $i + 1)
nam=`basename $line`
sbatch -A $proj --mail-user $mail -e $out/$line.err -o $out/$line.out -J Trim-$line -d afterok:$(expr $id + $i)\
 $UPSCb/pipeline/runTrimmomatic.sh $in/${line}_1.fq.gz $in/${line}_2.fq.gz $out
done

