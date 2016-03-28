#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
## $ -q nquintor.q
# $ -q miranda.q

# Run a list of files

#	set Geometry 		= 'Shell_VConst'
#	set Geometry 		= 'Wind2'
#	set Geometry 		= 'Static2_Wind2'
	set UseTburst		= 0

#	set Geometry 		= 'Static_Wind2'
#	set Geometry		= 'ThinShell'

#	set TagUV 			= '_UV_shield'
	set TagUV			= '_UV_none'
#	set TagUV 			= '_UV_shield'
	set TagLim 			= ''
	set EjectMode 		= 'all'
#	set Model 			= 'Bau05.zre10'
	set Method = 3
	set OutRad = '20.'
	set RootDir = '/home/aaorsi/LyaRT/'

#	foreach Model ('Bau05.zre10')
#	

	set MaxSlot	 = 6		# The actual number of running jobs will be MaxSlots * Num.Models * Num. Redshift, so be careful

	foreach Geometry ('Shell_VConst' 'ThinShell')
		if ($Geometry == 'Shell_VConst') then
			set TagUV	= '_UV_shield'
			set MFact		= '1.00'
			set RFactdisk 	= '0.500'
			set RFactbulge	= '0.500'
			set VFact		= '1.00'
		else
			set TagUV	= '_UV_none'
			set MFact		= '0.10'
			set RFactdisk 	= '2.000'
			set RFactbulge	= '2.000'
			set VFact		= '2.00'
		endif
	
		foreach Model (Bau05.zre10 Bau05.zre10_alphaesc0.0_resc0.2 Bau05.zre10_alphaesc0.0_resc0.5 Bau05.zre10_alphaesc1.0_resc0.1 Bau05.zre10_alphaesc1.0_resc0.2 Bau05.zre10_alphaesc1.5_resc0.2)
#		foreach Model (Bau05.zre10)

			foreach redshift (0.2 3.0 6.6)


				if ($Geometry == 'ThinShell' || $Geometry == 'Static_Wind2' || $Geometry == 'Static2_Wind2' || $Geometry == 'Static_Shell2' ) then
					set model = $Geometry$TagLim'_Method'$Method$TagUV'_mfact_'$MFact'_rfactdisk'$RFactdisk'_rfactbulge'$RFactbulge'_vfact'$VFact'_'$Model'_z'$redshift
				else
					set model = $Geometry$TagLim'_Method'$Method$TagUV'_Mdot_'$EjectMode'_rfactdisk'$RFactdisk'_rfactbulge'$RFactbulge'_vfact'$VFact'_outrad'$OutRad'_'$Model'_z'$redshift
	
				endif
				

				set Dir = $RootDir'/data/Params/'$model'/'
				set OutDir = $RootDir'/out/short/'$model'/'
				mkdir -p $OutDir
				set FName = $Dir'file_list'
	
				set nruns = `wc -l $FName`
		
				qsub -t 1-$nruns[1] -tc $MaxSlot -N 'LyaRT_'$Geometry -j y -o 'logfiles/' runLyaRT_list.csh $Dir $FName

			end
		end
	end
