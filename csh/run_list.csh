#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
# $ -q nquintor.q

# Run a list of files

#foreach mode (expanding static)
#	set Geometry 		= 'HomSphere'
	set Geometry 		= 'Shell_VConst'
#	set Geometry 		= 'HomShell'
#	set Geometry = 'HomSphere'
	set TagUV 			= '_UV_shield'
#	set TagUV 			= '_UV_none'
#	set TagUV 			= '_incUV'
#	set TagLim = '_incLim'
	set TagLim 			= '_noLim'
	set EjectMode 		= 'all'
	set Model 			= 'Bau05.zre10'
#	set Model = 'Bow06.zre10'
#	set Model = 'newmodel'

	
#	set TagUV = '_noUV'	
#set Geometry = 'HomSphere'

	set Method = 3
	set OutRad = '20.'

	foreach redshift (0.2 3.0 5.7 6.6)
#	foreach redshift (6.6)

#	set redshift = '0.3'
	set MFact = '1.00'	# Irrelevant if running the SN Outflows model
	set VFact = '0.50'
	set RFact = '1.00'

#	if ($Geometry == 'Shell_VConst') then
#		set model = $Geometry'_Method'$Method$TagUV'_'$mode'SN_Outflow_rfact'$RFact'_vfact'$VFact'_outrad'$OutRad'_Bau05.zre10_z'$redshift
		set model = $Geometry$TagLim'_Method'$Method$TagUV'_Mdot_'$EjectMode'_rfact'$RFact'_vfact'$VFact'_outrad'$OutRad'_'$Model'_z'$redshift
	if ($Geometry == 'HomShell') then
		set model = $Geometry$TagLim'_Method'$Method$TagUV'_mfact'$MFact'_rfact'$RFact'_vfact'$VFact'_'$Model'_z'$redshift
	endif

	set Dir = '/data/rw16/aaorsi/LyaRT/data/Params/'$model'/'
	set OutDir = '/data/rw16/aaorsi/LyaRT/out/short/'$model'/'
	mkdir -p $OutDir
	set FName = $Dir'file_list'

	set nruns = `wc -l $FName`
	
	qsub -t 1-$nruns[1] -N 'LyaRT_'$Geometry -j y -o 'logfiles/' runLyaRT_list.csh $Dir $FName

end

