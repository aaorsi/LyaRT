#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd


set ID = $1

set Seed0 = -20119
set isamp = $2
set tau0 = $3
set xcrit = $4
set FName = $5

@ Seed = $Seed0 + ($isamp * $isamp) * ($xcrit + 1)
echo $Seed
set ParFile = '../data/Params/'$ID'/'$FName$xcrit'.'$isamp'.par'
#echo $ParFile

LyaRT $ParFile $Seed


