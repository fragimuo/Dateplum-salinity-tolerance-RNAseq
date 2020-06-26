#!/bin/bash -l

set -ex

proj=u2018012
mail="francisco.gil.munoz@slu.se"
out=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/expressionQ
in=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/expressionQ
trans=/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity

sbatch -A $proj --mail-user $mail -e $out/ExpressionMatrix.err -o $out/ExpressionMatrix.out \
-J EMatrix $UPSCb/projects/persimmon/pipeline/runExpressionMatrixManual.sh $in $trans 
  
# Every variable passed after the .sh file is sended to the other script. They must be catched in the run script
# with the written order as $1, $2,... $n

