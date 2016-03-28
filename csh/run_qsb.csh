#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
## $ -q nquintor.q
# $ -q miranda.q

# Run a list of files

	set OutRad = '20.'
	set RootDir = '/home/aaorsi/LyaRT/'

#	foreach Model ('Bau05.zre10')
#	

	set MaxSlot	 = 24		# The actual number of running jobs will be MaxSlots * Num.Models * Num. Redshift, so be careful

	foreach Geometry ('Shell_VConst' 'ThinShell')
		if ($Geometry == 'Shell_VConst') then

			set mfq	 = '1.000'
			set vfq  = '1.000'
			set rfq  = '0.015'
			set mfsb = '1.000'
			set vfsb = '1.000'
			 
			set p0sb = '0.007'
			set p1sb = '2.700'

			set TagUV = '_UV_shield'

		else 
			set mfq	 = '0.100'
			set vfq  = '1.000'
			set rfq  = '0.400'
			set mfsb = '0.100'
			set vfsb = '1.000'
				 
			set p0sb = '0.200'
			set p1sb = '1.050'
	
			set TagUV = '_UV_none'
		endif
	
		foreach Model (Bau05.zre10.vcut30) # Bau05.zre10_alphaesc0.0_resc0.2 Bau05.zre10_alphaesc0.0_resc0.5 Bau05.zre10_alphaesc1.0_resc0.1 Bau05.zre10_alphaesc1.0_resc0.2 Bau05.zre10_alphaesc1.5_resc0.2)
#		foreach Model (Bau05.zre10)

			foreach redshift (0.2 3.0 4.5 5.7 6.6)


				set model = $Geometry$TagUV'_qmf'$qmf'_qrfd'$qrfd'_qrfb'$qrfb'_qvf'$qvf'_sbmf'$sbmf'_sbrfd'$sbrfd'_sbrfb'$sbrfb'_sbvf'$sbvf'_'$Model'_z'$redshift				

				set Dir = $RootDir'/data/Params/'$model'/'
				set OutDir = $RootDir'/out/short/'$model'/'
				mkdir -p $OutDir
				set FName = $Dir'file_list'
	
				set nruns = `wc -l $FName`
		
				qsub -t 1-$nruns[1] -tc $MaxSlot -N 'LyaRT_'$Geometry -j y -o 'logfiles/' runLyaRT_list.csh $Dir $FName

			end
		end
	end
