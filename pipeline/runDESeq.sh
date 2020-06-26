#! /bin/bash -l
#SBATCH -p core -n 1 -t 2-00:00:00
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
#


set -e

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
/usr/local/bin/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix $1/Trinity_trans.isoform.counts.matrix \
--samples_file $2 \
--method DESeq2 \
--output $3