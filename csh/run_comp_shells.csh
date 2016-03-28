#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
# $ -q mir_16.q


#	set model = 'comp_shells'
	set model = 'OutflowGrid'

	set LargeGrid = 1

	set MaxSlot = 100

#	set NoScattering = '_NoScatter'
	set NoScattering = ''

	set nh    = '22.0'
	set xcrit = '3'

	set Dir = '/home/aaorsi/LyaRT/data/Params/'$model'/'
#	set OutDir = '/home/aaorsi/LyaRT/out/short/'$model'/'
#	mkdir -p $OutDir

	if ($LargeGrid == 1) then
#		set FName = $Dir'file_list_large_xcrit'$xcrit$NoScattering
		set FName = $Dir'extra_list'
	else
		set FName = $Dir'file_list_nh'$nh'_xcrit'$xcrit
	endif
	
	set nruns = `wc -l $FName`
	
	qsub -t 1-$nruns[1] -tc $MaxSlot -N 'Shell_comp' -j y -o 'logfiles/' runLyaRT_list.csh $Dir $FName


