#$ -S /bin/tcsh -f
#PBS -o $PBS_O_WORKDIR/logfiles
#PBS -e $PBS_O_WORKDIR/logfiles

module purge
module load gsl-1.16

set Dir = $ARG1
set Seed0 = -69201
set id = $PBS_ARRAYID
set FName = $ARG2
set OutDir = $ARG3
set list = `cat $FName`

set ParFile = $Dir$list[$id]
echo $ParFile

set OFile = `echo $list[$id] | awk '{print substr($0,0,15)}'` 

set Outfile = $OutDir/$OFile

if (-f $Outfile) then
	echo $Outfile exists	
else
	@ Seed = $Seed0 + $id 
	echo $Seed
	echo $PBS_O_WORKDIR
	cd $PBS_O_WORKDIR  
	
	echo $Outfile does not exists
	/home/aaorsi/LyaRT/src/LyaRT $ParFile $Seed

endif



