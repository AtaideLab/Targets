#!/bin/bash

# Obtain list of IS from filtered BLAST output
filtered=Output/IS110_complete_Ident95_E0.bacterial_archaeal.blast.filtered
list=${filtered%.*}.list
awk 'NR>1 {print $1}' ${filtered} | sort | uniq > $list

# For each IS in list, concatenate all flanks into multi-flank-fasta
# and cleanup the temp fasta
indir=Output/EVEN-MORE-PARALLEL_Flanking_fastas_Ident95_E0
outdir=${indir}/400bp_concatentated_multifastas
mkdir -p $outdir

while read IS
do
	printf "Concatenating flanks for ${IS}\n"
	out=${outdir}/${IS}_400bp_flanks.fasta
	
	cat ${indir}/${IS}/${IS}*temp > ${out}
	
	rm -rf ${indir}/${IS}/${IS}*temp

done < $list

