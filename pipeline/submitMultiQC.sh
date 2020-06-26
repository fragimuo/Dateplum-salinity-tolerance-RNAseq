#!/bin/bash -l

set -e
# use set -ex for debugging instead

proj=u2018012
mail=francisco.gil.munoz@slu.se
in=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/multiqc

if [ ! -d $out ]; then
	mkdir -p $out
fi

module load bioinfo-tools multiqc

sbatch --mail-user=$mail -o $in/multiqc.out -e $in/multiqc.err \
-A $proj $UPSCb/pipeline/runMultiQC.sh $in $out
