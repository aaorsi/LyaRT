#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
## $ -q nquintor.q
## $ -q mir_16.q

set Dir = $1
set Seed0 = -11959
set id = $SGE_TASK_ID
set FName = $2

set list = `cat $FName`

set ParFile = $Dir$list[$id]
echo $ParFile

@ Seed = $Seed0 + $id 
echo $Seed

../src/LyaRT $ParFile $Seed


