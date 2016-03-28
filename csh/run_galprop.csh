#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
#$ -pe quintor 1


#Run a list of files

	set model = 'varying_galprops_incUV'
	set Dir = '/data/rw16/aaorsi/LyaRT/data/Params/'$model'/'

	set FName = $Dir'file_list'

	@ i = 0	
	foreach line (`cat $FName`)
				qsub -N 'LyaRT' -e 'logfiles/' -o 'logfiles/' runLyaRT_list.csh $i $Dir$line
			endif
		@ i++
	end


	
