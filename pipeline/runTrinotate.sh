#! /bin/bash -l
#SBATCH -p core -n 1 -t 2-00:00:00
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
#
# module load bioinfo-tools
# module load trinity

set -e

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-trinotate-berlin2018.simg \
/usr/local/src/TransDecoder/TransDecoder.LongOrfs -t $1

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-trinotate-berlin2018.simg \
/usr/local/src/TransDecoder/TransDecoder.Predict -t $1