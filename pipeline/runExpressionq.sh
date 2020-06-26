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
	Usage: $0 [options] <forward file> <reverse file> <transcriptome fasta> <temporal directory> <output directory>
  
  Options: 
    -p prepare the reference (default not to)
"
# load the helper (has functions that summarizes some functions)
source $UPSCb/src/bash/functions.sh

# default

PREP=
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

if [ "$#" -ne "5" ]; then
  abort "Number of variables different than expected"
fi

if [ ! -f $1 ]; then
  abort "Forward file not found or not a fastq file"
fi

if [ ! -f $2 ]; then
  abort "Reverse file not found or not a fastq file"
fi

if [ ! -f $3 ]; then
  abort "Transcriptome file not found or not a fasta file"
fi

if [ ! -d $4 ]; then
  abort "Temporal folder not found or not available"
fi

if [ ! -d $5 ]; then
  abort "Output folder not found or not available"
fi



#-d is a directory

#Unzipping files
fwd=$4/$(basename ${1/.gz/})
rev=$4/$(basename ${2/.gz/})
gunzip -c $1 > $fwd 
gunzip -c $2 > $rev

# run
singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
/usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl \
--seqType fq --left $fwd --right $rev --transcripts $3 \
--est_method salmon --trinity_mode --output_dir $5/$(basename ${1//_1.fq.gz/}) --thread_count $CPU $PREP

rm $fwd $rev

