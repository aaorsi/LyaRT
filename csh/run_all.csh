#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd


set Geom = 'HomSphere'
set FName = 'nc_1000_vel_0_xcrit'

set NSamp = (0)
@ it = 1
foreach tau0 (4.00)
	foreach xcrit (0)
		set ID = $Geom'_Tau'$tau0

		@ i = 0
		while ($i < $NSamp[$it])
			qsub -N $ID -e '../OutCDM/' -o '../OutCDM/' runLyaRT.csh $ID $i $tau0 $xcrit $FName
#			runLyaRT.csh $ID $i $tau0 $xcrit $FName
			@ i++
		end
	end
	@ it++
end



set Geom = 'HomSphere'

./../out/short/rm_out.csh
set tau0 = 7.06
@ it = 1
#set NSamp = (10 20 15 5)
set ncells = 1000
set NSamp = (15 15 10 10)
set xcrit = (3 3 3 3)

foreach ncells (1000)

	@ it = 1
	foreach vel (0 20 200 2000)

		set FName = 'nc_'$ncells'_vel_'$vel'_xcrit'
		set ID = $Geom'_Tau'$tau0

		@ i = 0
		while ($i < $NSamp[$it])
			qsub -N LyaRT -e '../OutCDM/' -o '../OutCDM/' runLyaRT.csh $ID $i $tau0 $xcrit[$it] $FName
			@ i++
		end
	    @ it++
	end
end



