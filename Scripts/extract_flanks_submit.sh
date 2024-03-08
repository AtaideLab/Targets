#!/bin/bash

#########################################################
#
# License: https://github.com/AtaideLab/Targets/blob/main/LICENSE
#
# Author/s: Cali Willet; cali.willet@sydney.edu.au
# Sydney Informatics Hub, The University of Sydney
#
#########################################################


#########################################################
# FILTERING PARAMETERS: CAN BE ADJUSTED BY USER OR LEAVE AS DEFAULT
#########################################################


# Path to your formatted BLAST non-redundant nucleotide database:
nt=/g/data/er01/NCBI/preformatted_2024-02-19/nt


# Prefix name of the multi-fasta  (omit .fasta suffix) 
#  containing all of the IS sequences you want to BLAST 
dataset=IS110_complete 


# Filter name:
# This is the name assigned within the filter blast step
# based on percent identity and E value thresholds set 
filter_name=Ident95_E0 


# Flank size to be extracted, ie take this many bp from 
# either side of the IS site from filtered BLAST hits
flank_size=200


#########################################################

sides=(left right)
script=Scripts/extract_flanks.pbs

for side in ${sides[@]}
do
	if [[ $1 == "test" ]]
	then
		# test/demo run:
		echo Running 2 IS demo: creating $side flanks
	
		nt=nt/IS_Targets_demo_db
		dataset=IS_2sequence_demo
	
		bash ${script} ${side} ${nt} ${dataset} ${filter_name} ${flank_size}
	
	else
		qsub -N ${side}-fl \
			-o ./PBS_logs/extract_flanks_${side}.o \
			-e ./PBS_logs/extract_flanks_${side}.e \
			-v side="${side}",nt="${nt}",dataset="${dataset}",filter_name="${filter_name}",flank_size="${flank_size}" \
			${script}
	fi
done
	

	

