#!/bin/bash -l

set -ex

proj=u2018012
tmp=/mnt/picea/tmp/fmunoz/expq
mail="francisco.gil.munoz@slu.se"
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/expressionQ
in=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/raw
trans=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity/Trinity.fasta

# Creating directories

if [ ! -d $tmp ]; then
	mkdir -p $tmp
fi

if [ ! -d $out ]; then
	mkdir -p $out
fi

# Loop for Expression quantification

i=1 #i=1 for indexing
for f in `find $in -name "*_[1,2].fq.gz" `; do echo "${f//_[1,2].fq.gz/}" \
; done | sort | uniq | while read line;
do 
  fnam=`basename $line`
  sbatch -A $proj --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out \
  -J smr-$fnam -d afterok:$(expr $id + $i)\
  $UPSCb/projects/persimmon/pipeline/runExpressionq.sh ${line}_1.fq.gz \
  ${line}_2.fq.gz $trans $tmp $out 
  
# Every variable passed after the .sh file is sended to the other script. They must be catched in the run script
# with the written order as $1, $2,... $n

done
