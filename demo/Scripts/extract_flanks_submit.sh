#!/bin/bash

#########################################################
#
# ADD LICENSE
#
# Author/s: Cali Willet; cali.willet@sydney.edu.au
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement: 
#       - See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST#acknowledgements
#
#
#########################################################

# simple wrapper script to submit 2 jobs with 1 variable (flank side):

sides=(left right)

# formatted BLAST non-redundant nucleotide database
nt=/g/data/er01/NCBI/preformatted_2024-02-19/nt

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
	

	

