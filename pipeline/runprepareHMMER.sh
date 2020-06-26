#! /bin/bash -l
#SBATCH -p core -n 20 -t 2-00:00:00
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
#
# module load bioinfo-tools
# module load trinity

set -e

module load bioinfo-tools
module load hmmer


hmmpress /mnt/picea/storage/reference/Pfam/Pfam31.0/Pfam-A.hmm 
          