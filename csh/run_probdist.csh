#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
## $ -q nquintor.q
# $ -q miranda.q

# Run a list of files

#	set Geometry 		= 'Shell_VConst'

	set Geometry	= 'HomSphere'
	set T			= '10000'	
	set xcrit		= '2.0'

	set model		= $Geometry'_T'$T'_xcrit'$xcrit'_fixedPosAndDir'

	set Dir = '/home/aaorsi/LyaRT/data/Params/ProbDist/'$model'/'
	set OutDir = '/home/aaorsi/LyaRT/out/short/ProbDist/'$model'/'
	mkdir -p $OutDir
	set FName = $Dir'file_list'

	set nruns = `wc -l $FName`
	echo $nruns	
	qsub -t 1-$nruns[1] -N 'LyaRT_'$Geometry -j y -o 'logfiles/' runLyaRT_list.csh $Dir $FName

