#! /bin/bash -l
#SBATCH -p core -n 48 --mem=300G -t 2-00:00:00
#SBATCH --mail-type=ALL --mail-user=francisco.gil.munoz@slu.se
#SBATCH -A u2018012
#
# module load bioinfo-tools
# module load trinity

set -ex

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
Trinity --seqType fq --max_memory 300G --samples_file $UPSCb/projects/persimmon/pipeline/samples.txt \
--CPU 30 --output /mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity \
--SS_lib_type RF 
