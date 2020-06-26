#!/bin/bash -l

## error and verbose
set -ex

# modules
module load bioinfo-tools hmmer

## check for the UPSCb env. var.
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
    usage
fi

## default args

proj=u2018012
mail="francisco.gil.munoz@slu.se"
in=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/hmmer/tmp/Trinity.fasta.transdecoder.pep
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/hmmer/hmmer
cpu=8

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch -A $proj --mail-user $mail -p all --array=1-400 -c $cpu -e $out/Trinity.prot.fasta.%a.err -o $out/Trinity.prot.fasta.%a.out \
$UPSCb/projects/persimmon/pipeline/runHMMER.sh $cpu $in $out 

## then
echo "Combine the results once done: cat Potra01-mRNA.blt.* > Potra01-mRNA.xml"
