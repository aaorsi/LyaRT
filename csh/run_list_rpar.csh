#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
## $ -q nquintor.q
# $ -q miranda.q

# Run a list of files

#	set Geometry 		= 'Shell_VConst'
#	set Geometry 		= 'Wind2'
#	set Geometry 		= 'Static2_Wind2'
	set UseTburst		= 1

	set Geometry 		= 'Static_Wind2'
#	set Geometry		= 'ThinShell'

#	set TagUV 			= '_UV_shield'
	set TagUV			= '_UV_none'
#	set TagUV 			= '_UV_shield'
	set TagLim 			= '_noLim'
	set EjectMode 		= 'all'
	set Model 			= 'Bau05.zre10'
	set Method = 3
	set OutRad = '20.'

	foreach Model ('Bau05.zre10')

#		foreach redshift (0.2 3.0 5.7 6.6)
	foreach redshift (0.2 3.0 6.6)

#	set redshift = '0.3'
			set MFact		= '1.00'
			set RFactdisk 	= '1.500'
			set RFactbulge	= '1.500'
			set VFact		= '1.00'

			if ($Geometry == 'ThinShell' || $Geometry == 'Static_Wind2' || $Geometry == 'Static2_Wind2' || $Geometry == 'Static_Shell2' ) then
				set model = $Geometry$TagLim'_Method'$Method$TagUV'_mfact_'$MFact'_rfactdisk'$RFactdisk'_rfactbulge'$RFactbulge'_vfact'$VFact'_'$Model'_z'$redshift
			else
				set model = $Geometry$TagLim'_Method'$Method$TagUV'_Mdot_'$EjectMode'_rfactdisk'$RFactdisk'_rfactbulge'$RFactbulge'_vfact'$VFact'_outrad'$OutRad'_'$Model'_z'$redshift

			endif
				
			if ($UseTburst == 1) then
				set model = $Geometry'_tburst'$TagLim'_Method'$Method$TagUV'_mfact_'$MFact'_vfact'$VFact'_'$Model'_z'$redshift



			set Dir = '/data/rw16/aaorsi/LyaRT/data/Params/'$model'/'
			set OutDir = '/data/rw16/aaorsi/LyaRT/out/short/'$model'/'
			mkdir -p $OutDir
			set FName = $Dir'file_list'

			set nruns = `wc -l $FName`
	
			qsub -t 1-$nruns[1] -N 'LyaRT_'$Geometry -j y -o 'logfiles/' runLyaRT_list.csh $Dir $FName

		end
	end

