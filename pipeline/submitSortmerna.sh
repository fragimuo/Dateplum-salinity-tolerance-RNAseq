#!/bin/bash -l

set -e

proj=u2018012
tmp=/mnt/picea/tmp/fmunoz
mail="francisco.gil.munoz@slu.se"
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/sortmerna
in=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/raw

module load bioinfo-tools
module load sortmerna

for f in `find $in -name "*_[1,2].fq.gz" `; do echo "${f//_[1,2].fq.gz/}" ; done | sort | uniq | while read line;
do 
fnam=`basename $line`
sbatch -A $proj --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out -J smr-$fnam $UPSCb/pipeline/runSortmerna.sh -m \
$out $tmp  ${line}_1.fq.gz ${line}_2.fq.gz
done

