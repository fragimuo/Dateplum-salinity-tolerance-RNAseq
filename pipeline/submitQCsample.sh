#!/bin/bash -l

set -ex

proj=u2018012
mail="francisco.gil.munoz@slu.se"
in=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity


sbatch -A $proj --mail-user $mail -e $in/QCsample.err -o $in/QCsample.out \
-J QC $UPSCb/projects/persimmon/pipeline/runQCsample.sh $in
  
# Every variable passed after the .sh file is sended to the other script. They must be catched in the run script
# with the written order as $1, $2,... $n

