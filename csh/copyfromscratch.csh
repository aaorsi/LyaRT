#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd

set NNodes = 64
set scratchdir = '/scratch/aaorsi/LyaRT/out/short/OutflowGrid/'
set OutDir 	   = '/fdg/aaorsi/LyaRT/out/short/OutflowGrid/'

@ i = 1

	while ($i <= $NNodes) 
		if ($i < 10) then 
			set inode = '0'$i
		else
			set inode = $i
		endif


		echo geryon$inode
		ssh 'geryon'$inode cp  $scratchdir/'*' $OutDir
		echo 'done!'
		


		@ i++

	end


