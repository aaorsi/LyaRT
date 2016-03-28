#$ -S /bin/tcsh -f

module load gsl-1.16

set ParDir 	 = '/home/aaorsi/LyaRT/data/Params/'
set ODir 	 = '/home/aaorsi/LyaRT/out/short/'
set GridType = 'MoR_V_Z'

set Size     = '0.1x0.1x0.1'

set	OutGridName = $GridType'.'$Size
#set Model = 'Bau05.zre10.vcut30'
set Model = 'Uber'

set MaxSlot	= 64
set PreName = ''

#foreach Geometry ('Shell_VConst')# 'ThinShell')
#foreach Geometry ('Shell_VConst')
foreach Geometry ('ThinShell')

	if ($Geometry == 'Shell_VConst') then

		set mfq	 = '1.000'
		set vfq  = '1.000'
		set rfq  = '0.250'
		set mfsb = '1.000'
		set vfsb = '1.000'
			 
		set p0sb = '0.082'
		set p1sb = '0.636'

		set TagUV = '_UV_shield'
		set rf = (0.100 0.600 0.800)	# one for each redshift, see below

	else 
		set mfq	 = '0.100'
		set vfq  = '1.000'
#		set rfq  = '1.200'
		set rfq  = '2.000'
		set mfsb = '0.100'
		set vfsb = '1.000'
			 
#		set p0sb = '0.148'
#		set p1sb = '1.086'

		set p0sb = '0.200'
		set p1sb = '1.000'

		set TagUV = '_UV_none'
		set rf = (0.300 2.000 1.800)
	endif

#	foreach redshift (0.2 3.0 4.5 5.7 6.6 7.3)
	@ i = 1
	foreach redshift (4.5 5.7)# 4.5 5.7)
		
#		if ($Model == 'Uber') then
#			set NameGrid = $PreName'qmf'$mfq'qvf'$vfq'sbmf'$mfsb'sbvf'$vfsb'rf'$rf[$i]'_'$Geometry'_'$Model'_z'$redshift
#		else
			set NameGrid = $PreName'qmf'$mfq'qvf'$vfq'qrf'$rfq'sbmf'$mfsb'sbvf'$vfsb'p0sb'$p0sb'p1sb'$p1sb'_'$Geometry'_'$Model'_z'$redshift
#		endif
		set GridDir	 = 'Grid/'$OutGridName'/'$NameGrid'/'

		set PDir	= $ParDir'/'$GridDir
		set OutDir	 = $ODir'/'$GridDir

		mkdir -p $OutDir

		set GridFile = $PDir'/grid_list'

		set nruns 	 = `wc -l $GridFile`

		echo $redshift $Geometry $nruns[1] 


#		original:
#		qsub -t 1-$nruns[1] -tc $MaxSlot -N 'LyaGrid_'$Size -j y -o 'logfiles/' runLyaRT_grid.csh $PDir $GridFile
#		qsub -t 1-$nruns[1]%10 -j oe -N LyaRT_z$redshift -v ARG1=$PDir,ARG2=$GridFile ./runLyaRT_grid.csh
		qsub -t 1-$nruns[1]%$MaxSlot -k oe -N LyaRT_z$redshift -v ARG1=$PDir,ARG2=$GridFile,ARG3=$OutDir ./runLyaRT_grid.csh
#		qsub -t 1-100 -k oe -N LyaRT_z$redshift -v ARG1=$PDir,ARG2=$GridFile ./runLyaRT_grid.csh

#		if ($redshift == 3.0) then
#			qsub -t 14400-$nruns[1] -tc $MaxSlot -N 'LyaGrid_'$Size -j y -o 'logfiles/' runLyaRT_grid.csh $PDir $GridFile
#		else if ($redshift == 4.5) then
#			qsub -t 12300-$nruns[1] -tc $MaxSlot -N 'LyaGrid_'$Size -j y -o 'logfiles/' runLyaRT_grid.csh $PDir $GridFile
#		else if ($redshift == 5.7) then
#			qsub -t 11500-$nruns[1] -tc $MaxSlot -N 'LyaGrid_'$Size -j y -o 'logfiles/' runLyaRT_grid.csh $PDir $GridFile
#		else if ($redshift == 7.3) then 		 
#			qsub -t 10000-$nruns[1] -tc $MaxSlot -N 'LyaGrid_'$Size -j y -o 'logfiles/' runLyaRT_grid.csh $PDir $GridFile
#		endif 
		@ i++
	end
end
