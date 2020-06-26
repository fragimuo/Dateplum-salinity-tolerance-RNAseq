#! /bin/bash -l
#SBATCH -p core -n 1 -t 2-00:00:00
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
#
# module load bioinfo-tools
# module load trinity

set -e


singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-trinotate-berlin2018.simg \
 signalp -f short -n $2/signalp.out $1.$SLURM_ARRAY_TASK_ID > $2/sigP.$SLURM_ARRAY_TASK_ID.log-t

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-trinotate-berlin2018.simg \
 /usr/local/bin/tmhmm --short < $1.$SLURM_ARRAY_TASK_ID > $2/tmhmm.$SLURM_ARRAY_TASK_ID.out