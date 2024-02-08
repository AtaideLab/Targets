#!/bin/bash

# Takes the split output from extract_flanks_parallel.pbs
# and combines them into a per-IS multi-fasta
# Number of outputs should be equal to the number of IS
# with >=1 hit passing the specified filter 

#########################################################
# PARAMETERS: TO BE ADJUSTED BY USER
#########################################################

# Output file aming used in previous steps: 
filterName=Ident95_E0
dataset=IS110_complete

#########################################################

indir=./Output/Flanking_fastas_${filterName}
report=./Output/${dataset}_${filterName}.bacterial.blast.report

awk 'NR>1' ${report} | while read LINE
do
	is=$(echo $LINE | awk '{print $1}')
	fam=$(echo $LINE | awk '{print $2}')
	group=$(echo $LINE | awk '{print $3}')
	
	if [[ $group =~ ^NA$ ]]
	then
		group=unknown
	fi
	
	hits=$(echo $LINE | awk '{print $5}')
	
	is_id=${is}_${fam}_${group}
	
	if [[ $hits -gt 0 ]]
	then
		files=$(ls ${indir}/${is_id}_400bp_flanks.fasta-*)
		cat $files > ${indir}/${is_id}_400bp_flanks.fasta
	fi
done	
	
