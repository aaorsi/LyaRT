#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
# $ -q nquintor.q

# Run a list of files


set TagUV	 = '_UV_none'	 
set ParDir 	 = '/data/rw16/aaorsi/LyaRT/data/Params/Grid/'
set ODir 	 = '/data/rw16/aaorsi/LyaRT/out/short/Grid/'
set GridType = 'MoR_V_Z'
set Size     = '0.1x0.1x0.1'

set vf		= '0.50'
set rf		= '1.00'

set PDir	= $ParDir$GridType'/'$GridType'.'$Size$TagUV'/'
set OutDir	 = $ODir$GridType'/'$GridType'.'$Size$TagUV'/'

mkdir -p $OutDir
set GridFile = $PDir'/grid_list.vf'$vf'_rf'$rf

set nruns 	 = `wc -l $GridFile`

qsub -t 1-$nruns[1] -N 'LyaGrid_'$Size -j y -o 'logfiles/' runLyaRT_grid.csh $PDir $GridFile


