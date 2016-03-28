#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd


set ID = $1

set Seed0 = 21255
set isamp = $2
set tau0 = $3
set xcrit = $4
set atauf = $5

set Seed = `echo $Seed0 $isamp $xcrit $atauf |awk '{print $1 + ($2 * $2) * ($3 + 1) - ($4 * 100)}'`
set ParFile = '../data/Params/'$ID'/afac'$atauf'xcrit'$xcrit'_'$isamp'.par'

LyaRT $ParFile $Seed

#echo $ParFile $Seed

