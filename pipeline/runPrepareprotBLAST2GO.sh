#! /bin/bash -l
#SBATCH -p core -n 1 -t 2-00:00:00
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
module load R

./prepareprotBLAST2GO.R