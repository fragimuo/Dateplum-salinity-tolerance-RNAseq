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
Usage: $0 <Transcripts matrix directory>

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

if [ "$#" -ne "1" ]; then
abort "Number of variables different than expected"
fi

if [ ! -d $1 ]; then
abort "Expression matrix not found"
fi


#-d is a directory


singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
/usr/local/bin/trinityrnaseq/Analysis/DifferentialExpression/PtR \
 -m $1/Trinity_trans.isoform.counts.matrix -s $1/samples.txt \
 --log2 --compare_replicates 
 
singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
/usr/local/bin/trinityrnaseq/Analysis/DifferentialExpression/PtR \
      -m $1/Trinity_trans.isoform.TMM.EXPR.matrix -s $1/samples.txt \
      --log2 --sample_cor_matrix
 
singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
/usr/local/bin/trinityrnaseq/Analysis/DifferentialExpression/PtR\
       -m $1/Trinity_trans.isoform.TMM.EXPR.matrix -s $1/samples.txt\
       --log2 --CPM --prin_comp 3