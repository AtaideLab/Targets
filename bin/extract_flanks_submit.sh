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

script=bin/extract_flanks.pbs

for side in ${sides[@]}
do
	qsub -N ${side}-fl \
		-o ./PBS_logs/extract_flanks_${side}.o \
		-e ./PBS_logs/extract_flanks_${side}.e \
		-v side="${side}" \
		${script}
done
	

	

