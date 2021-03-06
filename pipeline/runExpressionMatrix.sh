#! /bin/bash -l
#SBATCH -p core -n 30 --mem=50G -t 2-00:00:00
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
#
# module load bioinfo-tools
# module load trinity

set -ex

#Calling for options 
# $@
# $1
# $2
# $3
# $4

# usage text

export USAGETXT=\
"
	Usage: $0  <expression quantification folder> <transcriptome folder> 
  
"
# load the helper (has functions that summarizes some functions)
source $UPSCb/src/bash/functions.sh

# default

CPU=30

# Options. We are defining the --prep_reference Trinity option as a variable to be written in the command

while getopts hp option
do
  case "$option" in
      h) usage;;
      p) PREP="--prep_reference";;
      ?) usage;;
  esac
done
shift `expr $OPTIND - 1`


# Defining variables and testing for them ($number is each variable call in the order
# of they were defined after calling the script)

if [ "$#" -ne "2" ]; then
  abort "Number of variables different than expected"
fi

if [ ! -d $1 ]; then
  abort "Quantification files not found"
fi

if [ ! -d $2 ]; then
  abort "Transcriptome folder not found"
fi


#-d is a directory



# run
singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
/usr/local/bin/trinityrnaseq/util/abundance_estimates_to_matrix.pl \
--est_method salmon  --out_prefix Trinity_trans \
--name_sample_by_basedir \
for f in `find $1 -name "*.out" `; do echo "${f//.out/}/quant.sf \ " \
; done | sort | uniq \
--gene_trans_map $2/Trinity.fasta.gene_trans_map
