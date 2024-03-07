#!/bin/bash

#########################################################
#
# License: https://github.com/AtaideLab/Targets/blob/main/LICENSE
#
# Author/s: Cali Willet; cali.willet@sydney.edu.au
# Sydney Informatics Hub, The University of Sydney
#
#########################################################

# simple wrapper script to submit 2 jobs with 1 variable (flank side):

sides=(left right)

# path to your formatted BLAST non-redundant nucleotide database:
nt=<filepath>/preformatted_2024-02-19/nt

# the prefix name of the multi-fasta  (omit .fasta/.fa suffix) containing all of the IS sequences you want to BLAST 
dataset=IS110_complete 

script=Scripts/extract_flanks.pbs

for side in ${sides[@]}
do
	if [[ $1 == "test" ]]
	then
		# test/demo run:
		echo Running 2 IS demo: creating $side flanks
	
		nt=nt/IS_Targets_demo_db
		dataset=IS_2sequence_demo
	
		bash ${script} ${side} ${nt} ${dataset}
	
	else
		qsub -N ${side}-fl \
			-o ./PBS_logs/extract_flanks_${side}.o \
			-e ./PBS_logs/extract_flanks_${side}.e \
			-v side="${side}" \
			${script}
	fi
done
	

	

