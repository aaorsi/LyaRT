#$ -S /bin/tcsh
# Run a list of files

set ParDir 	 = '/home/aaorsi/LyaRT/'
set ODir 	 = '/home/aaorsi/LyaRT/'

#	set Geometry 		= 'Shell_VConst'
#	set Geometry 		= 'Wind2'
#	set Geometry 		= 'Static2_Wind2'
#	set PreName			= 'Feldmeier12'		# Leave blank otherwise
	set PreName			= 'Uber'		# Leave blank otherwise
#	set PreName			= 'Uber_powlaw'		# Leave blank otherwise
	set UseTburst		= 0

#	set Geometry 		= 'Static_Wind2'
#	set Geometry		= 'ThinShell'

#	set TagUV 			= '_UV_shield'
#	set TagUV			= '_UV_none'
#	set TagUV 			= '_UV_shield'
	set TagLim 			= ''
	set EjectMode 		= 'all'
#	set Model 			= 'Bau05.zre10'
#	set RootDir = '/VVV_data2/aaorsi/LyaRT/'

#	foreach Model ('Bau05.zre10')
#	

	set MaxSlot	 = 33		# The actual number of running jobs will be MaxSlots * Num.Models * Num. Redshift, so be careful

#	echo 'Sleeping for two long hours...'	
#	sleep 7200		# Sleep for two hours
	
	foreach Geometry ('ThinShell')
#	foreach Geometry ('Shell_VConst')
#	foreach Geometry ('Wind_VelProfile')
		if ($Geometry == 'Shell_VConst') then
	
		set mfq	 = '1.000'
		set vfq  = '1.000'
		set rfq  = '2.000'
		set mfsb = '1.000'
		set vfsb = '1.000'
			 
#		set p0sb = '0.010'
#		set p0sb = '0.014'
#		set p1sb = '2.152'
		
		set p0sb = '0.263'
		set p1sb = '0.921'

		set TagUV = '_UV_shield'
		
#		set rf = (0.030 0.700 1.000)	# one for each redshift, see below
#		set rf = (2.000 3.000)	# one for each redshift, see below

		else if ($Geometry == 'Wind_VelProfile') then
	
		set mfq	 = '1.000'
		set vfq  = '1.000'
		set rfq  = '2.000'
		set mfsb = '1.000'
		set vfsb = '1.000'
			 
#		set p0sb = '0.010'
		set p0sb = '0.200'
		set p1sb = '1.000'

		set TagUV = '_UV_none'

		else if ($Geometry == 'ThinShell')  then 

		set mfq	 = '0.100'
		set vfq  = '1.000'
		set rfq  = '2.000'
		set mfsb = '0.100'
		set vfsb = '1.000'
			 
#		set p0sb = '0.150'
#		set p0sb = '0.223'
#		set p1sb = '0.925'

		set p0sb = '0.200'
		set p1sb = '1.000'

		set TagUV = '_UV_none'
		
#		set rf = (0.100 1.000 1.300)

		endif
	
#		foreach Model (Bau05.zre10.vcut30) # Bau05.zre10_alphaesc0.0_resc0.2 Bau05.zre10_alphaesc0.0_resc0.5 Bau05.zre10_alphaesc1.0_resc0.1 Bau05.zre10_alphaesc1.0_resc0.2 Bau05.zre10_alphaesc1.5_resc0.2)
		foreach Model (Uber) 

# Bau05.zre10_alphaesc0.0_resc0.2 Bau05.zre10_alphaesc0.0_resc0.5 Bau05.zre10_alphaesc1.0_resc0.1 Bau05.zre10_alphaesc1.0_resc0.2 Bau05.zre10_alphaesc1.5_resc0.2)
#		foreach Model (Bau05.zre10)

			@ i = 1
			foreach redshift (2.2 3.0 5.7)
#			foreach redshift (3.0 5.7)

#				set model = $Geometry$TagUV'_qmf'$mfq'_qrfd'$rfq'_qrfb'$qrfb'_qvf'$qvf'_sbmf'$sbmf'_p0'$p0'_p1'$p1'_sbvf'$sbvf'_'$Model'_z'$redshift				
				if ($Model == 'Uber_Set1') then
					set model = $PreName$Geometry$TagUV'qmf'$mfq'qvf'$vfq'sbmf'$mfsb'sbvf'$vfsb'rf'$rf[$i]'_Uber_z'$redshift
				else 
				set model = $PreName$Geometry$TagUV'qmf'$mfq'qvf'$vfq'qrf'$rfq'sbmf'$mfsb'sbvf'$vfsb'p0sb'$p0sb'p1sb'$p1sb'_'$Model'_z'$redshift
				endif

				set Dir = $ParDir'/data/Params/'$model'/'
				set OutDir = $ODir'/out/short/'$model'/'
				mkdir -p $OutDir
				set FName = $Dir'file_list'
	
				set nruns = `wc -l $FName`
		
				qsub -t 1-$nruns[1]%$MaxSlot -k oe -N 'LyaRT_'$Geometry -v ARG1=$Dir,ARG2=$FName ./runLyaRT_list.csh 
#				echo 'qsub -t 1-'$nruns[1]' -tc '$MaxSlot '-N LyaRT_'$Geometry '-j y -o logfiles/ runLyaRT_list.csh '$Dir $FName
			
				@ i++

			end
		end
	end
