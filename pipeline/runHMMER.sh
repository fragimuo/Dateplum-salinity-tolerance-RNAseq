#! /bin/bash -l
#SBATCH -p all
#SBATCH -n 1
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
#
# module load bioinfo-tools
# module load trinity

set -e

module load bioinfo-tools
module load hmmer

fnam=$(basename ${2//.f*a*/})

hmmscan --cpu $1 --domtblout $3/TrinotatePFAM.$SLURM_ARRAY_TASK_ID.out \
          /mnt/picea/storage/reference/Pfam/Pfam31.0/Pfam-A.hmm \
          $2.$SLURM_ARRAY_TASK_ID > $3/$fnam.$SLURM_ARRAY_TASK_ID.pfam.log