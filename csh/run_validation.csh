#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd
# # $ -q nquintor.q
# # $ -q mir_16.q

# Run a list of files, different tests

set RemoveOldData	= 'no'


set Test_HomSlab 	= 1
set Test_HomSphere 	= 0
set Test_NScat		= 0
set Test_fesc		= 0
set Test_ExpSphere	= 0
set Test_ThinShell	= 0

set ParDir 	 = '/home/aaorsi/LyaRT/data/Params/validation/'
set OutDir 	 = '/home/aaorsi/LyaRT/out/validation/'

set MaxSlot	 = 100


if ($RemoveOldData == 'yes') then
	echo 'Removing everything from '$OutDir
	rm -r $OutDir
endif

if ($Test_HomSlab == 1) then

	set Geom	= 'HomSlab'
	set T 		= '10.00'
	foreach xcrit ('0.0' '3.0') 
		set NameDir = $Geom'_T'$T'_xcrit'$xcrit'/'	
		set ODir	= $OutDir$NameDir
		mkdir -p $ODir
		
		set Dir		= $ParDir$NameDir
	
		set ListFile = $Dir'/file_list'
	
		set nruns 	 = `wc -l $ListFile`

		echo $ODir
	
		qsub -t 1-$nruns[1] -tc $MaxSlot -N 'Lya_test_homslab' -j y -o 'logfiles/' runLyaRT_tests.csh $Dir $ListFile
	end
	
endif

if ($Test_HomSphere == 1) then

	set Geom	= 'HomSphere'
	set T 		= '10.00'
	set xcrit 	= '3.0' 
	set NameDir = $Geom'_T'$T'_xcrit'$xcrit'/'	
	set ODir	= $OutDir$NameDir
	mkdir -p $ODir
	set Dir		= $ParDir$NameDir

	set ListFile = $Dir'/file_list'

	set nruns 	 = `wc -l $ListFile`

	qsub -t 1-$nruns[1] -tc $MaxSlot -N 'Lya_test_homsphere' -j y -o 'logfiles/' runLyaRT_tests.csh $Dir $ListFile

endif

if ($Test_fesc == 1) then

	set Geom	= 'HomSlab'
	set T 		= '10.00'
	set xcrit 	= '3.0' 
	set NameDir = 'fesc_'$Geom'_T'$T'_xcrit'$xcrit'/'	
	set ODir	= $OutDir$NameDir
	mkdir -p $ODir
	set Dir		= $ParDir$NameDir

	set ListFile = $Dir'/file_list'

	set nruns 	 = `wc -l $ListFile`

	qsub -t 1-$nruns[1] -tc $MaxSlot -N 'Lya_test_homslab' -j y -o 'logfiles/' runLyaRT_tests.csh $Dir $ListFile

endif


if ($Test_ExpSphere == 1) then

	set Geom	= 'HomSphere'
	set T 		= '10000'
	set xcrit 	= '3.0' 
	set NameDir = 'Expanding_'$Geom'_T'$T'_xcrit'$xcrit'/'	
	set ODir	= $OutDir$NameDir
	mkdir -p $ODir
	set Dir		= $ParDir$NameDir

	set ListFile = $Dir'/file_list'

	set nruns 	 = `wc -l $ListFile`

	qsub -t 1-$nruns[1] -tc $MaxSlot -N 'Lya_test_exp' -j y -o 'logfiles/' runLyaRT_tests.csh $Dir $ListFile

endif

if ($Test_ThinShell == 1) then

	set Geom	= 'ThinShell'
	set T 		= '96897'
#	set T 		= '10000'
	set xcrit 	= '0.0' 
	set NameDir = $Geom'_T'$T'_xcrit'$xcrit'/'	
	set ODir	= $OutDir$NameDir
	mkdir -p $ODir
	set Dir		= $ParDir$NameDir

	set ListFile = $Dir'/file_list'

	set nruns 	 = `wc -l $ListFile`

	qsub -t 1-$nruns[1] -tc $MaxSlot -N 'Lya_thinshell' -j y -o 'logfiles/' runLyaRT_tests.csh $Dir $ListFile

endif




if ($Test_NScat == 1) then

	set Geom	= 'HomSlab'
	set T 		= '10.00'
	set xcrit 	= '0.0' 
	set NameDir = 'NScat'$Geom'_T'$T'_xcrit'$xcrit'/'	
	set ODir	= $OutDir$NameDir
	mkdir -p $ODir
	set Dir		= $ParDir$NameDir

	set ListFile = $Dir'/file_list'

	set nruns 	 = `wc -l $ListFile`

	qsub -t 1-$nruns[1] -tc $MaxSlot -N 'Lya_test_homslab' -j y -o 'logfiles/' runLyaRT_tests.csh $Dir $ListFile

endif


