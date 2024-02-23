#!/bin/bash

# Obtain list of IS from filtered BLAST output
filtered=Output/IS110_complete_Ident95_E0.bacterial_archaeal.blast.filtered
list=${filtered%.*}.list
awk 'NR>1 {print $1}' ${filtered} | sort | uniq > $list

list=list
# For each IS in list, concatenate all flanks into multi-flank-fasta
# and cleanup the temp dirs

dir=2seq_test

while read IS
do
	printf "Concatenating flanks for ${IS}\n"
	out=${dir}/${IS}_400bp_flanks.fasta
	cat ${dir}/${IS}/* > ${out}
done < $list

