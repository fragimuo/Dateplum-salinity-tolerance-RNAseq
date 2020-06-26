#! /bin/bash -l
#SBATCH -p core -n 1 -t 01:00:00
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
#
# module load bioinfo-tools
# module load trinity

set -ex

export USAGETXT=\
"
Usage: $0 <Transcripts matrix> <transcriptome fasta> <output directory>

"
# load the helper (has functions that summarizes some functions)
source $UPSCb/src/bash/functions.sh


# Options. We are defining the --prep_reference Trinity option as a variable to be written in the command

while getopts hp option
do
case "$option" in
h) usage;;
?) usage;;
esac
done
shift `expr $OPTIND - 1`


# Defining variables and testing for them ($number is each variable call in the order
# of they were defined after calling the script)

if [ "$#" -ne "3" ]; then
abort "Number of variables different than expected"
fi

if [ ! -f $1 ]; then
abort "Expression matrix not found"
fi

if [ ! -f $2 ]; then
abort "Transcriptome file not found or not a fasta file"
fi

if [ ! -d $3 ]; then
abort "Output folder not found or not available"
fi


#-d is a directory



singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
/usr/local/bin/trinityrnaseq/util/misc/contig_ExN50_statistic.pl \
$1 $2 | tee $3\ExN50.stats