#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
# #$ -q nquintor.q
# #$ -q mir_16.q


set UseScratch = 1
set scratchdir = '/scratch/aaorsi/LyaRT/out/short/OutflowGrid/'
if ($UseScratch == 1) then 
	mkdir -p $scratchdir
endif

set Dir = $1
set Seed0 = 7129101
set id = $SGE_TASK_ID
set FName = $2

set list = `cat $FName`

set ParFile = $Dir$list[$id]
echo $ParFile

@ Seed = $Seed0 + $id 
echo $Seed

../src/LyaRT $ParFile $Seed

