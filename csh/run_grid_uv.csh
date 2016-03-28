#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
# $ -q nquintor.q

# Run a list of files


set ParDir 	 = '/data/rw16/aaorsi/LyaRT/data/Params/Grid/'
set ODir 	 = '/data/rw16/aaorsi/LyaRT/out/short/Grid/'
set GridType = 'MoR_V_Z'

#set Model	 = 'ThinShell'
#set Model	 = 'Shell_VConst'
set Model	 = 'Wind2'

set Size     = '0.1x0.1x0.1'

set mf		= '1.00'
set vf		= '1.00'
set rfdisk	= '0.02'
set rfbulge	= '0.02'

#set TagUV   = '_UV_none'
set TagUV   = '_UV_shield'

#foreach z (0.2 3.0 5.7 6.6)
foreach z (0.2)

	set PDir	= $ParDir$Model'/'$GridType'/'$GridType'.'$Size$TagUV'_z'$z'/'
	set OutDir	 = $ODir$Model'/'$GridType'/'$GridType'.'$Size$TagUV'_z'$z'/'

	mkdir -p $OutDir
	set GridFile = $PDir'/grid_list.mf'$mf'vf'$vf'rfactdisk'$rfdisk'rfactbulge'$rfbulge

	set nruns 	 = `wc -l $GridFile`

	qsub -t 5500-$nruns[1] -N 'LyaGrid_'$Size -j y -o 'logfiles/' runLyaRT_grid.csh $PDir $GridFile
end
