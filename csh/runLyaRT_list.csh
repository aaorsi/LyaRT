#$ -S /bin/tcsh -f
#PBS -o $PBS_O_WORKDIR/logfiles
#PBS -e $PBS_O_WORKDIR/logfiles

module purge
module load gsl-1.16

set UseScratch = 0
set scratchdir = '/scratch/aaorsi/LyaRT/out/short/OutflowGrid/'
if ($UseScratch == 1) then 
	mkdir -p $scratchdir
endif

set Dir = $ARG1
set Seed0 = 7129101
set id = $PBS_ARRAYID
set FName = $ARG2

set list = `cat $FName`

set ParFile = $Dir$list[$id]
echo $ParFile

@ Seed = $Seed0 + $id 
echo $Seed

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR  
/home/aaorsi/LyaRT/src/LyaRT $ParFile $Seed

