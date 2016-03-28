#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd

set NSamp = 25

foreach tau0 (6.00)
	set ID = 'HomSlab_Tau'$tau0'_dust'
	foreach xcrit (3)
		foreach atauf (1.25 1.58 1.99 2.51 3.16 3.98 5.01 6.30 7.94 12.58)
			@ i = 0
			while ($i < $NSamp)
				qsub -N $ID -e '../OutCDM/' -o '../OutCDM/' runLyaRT_dust.csh $ID $i $tau0 $xcrit $atauf
#				runLyaRT_dust.csh $ID $i $tau0 $xcrit $atauf
				@ i++
			end
		end
	end
end



