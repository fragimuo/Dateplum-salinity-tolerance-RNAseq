#! /bin/bash -l
#SBATCH -p core -n 48 --mem=300G -t 2-00:00:00
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
#
# module load bioinfo-tools
# module load trinity

set -ex

#Calling for options 
#$@
# $1
# $2
# $3
# $4

# tests
if [ "$#" -ne "4" ]; 
  echo "Number of variables different than expected"
  exit 1
fi

if [ ! -f $1 ];
  echo
  exit 1
fi

if [ ! -f $2 ];
  echo
  exit 1
fi

if [ ! -d $3 ];
  echo
  exit 1
fi

if [ ! -d $4 ];
  echo
  exit 1
fi

#-d is a directory

#Unzipping files
fwd=$tmp/$(basename ${1/.gz/})
rev=$tmp/$(basename ${2/.gz/})
gunzip -c $1 > $fwd 
gunzip -c $2 > $rev

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
$TRINITY_HOME/util/align_and_estimate_abundance.pl \
        --seqType fq \
        --left ${line}_1.fq \ 
        --right ${line}_2.fq \
        --transcripts /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity/Trinity.fasta  \
        --est_method salmon  --aln_method bowtie \
        --trinity_mode \
        --output_dir $out \

        
        