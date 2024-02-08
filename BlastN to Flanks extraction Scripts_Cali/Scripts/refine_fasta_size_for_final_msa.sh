#!/bin/bash

#########################################################
#
# Platform: USyd Artemis HPC
#
# Description: 
#	- Create new fasta for each IS restricted to specified range
#	- Range for IS1111 = 1000 bp onwards
#	- Range for IS110 = 1-500 bp 
#
# Usage:
#       - bash Scripts/refine_fasta_size_for_final_msa.sh
#
# Compute resources:
#       - a few minutes on login node 
#
# Output:
#       - Input/<IS>_<BP>.fasta - fasta for each IS restricted to specified range
#
# Author/s: Cali Willet; cali.willet@sydney.edu.au
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement: 
#       - See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST#acknowledgements
#
#########################################################

#ISPsy7_IS110_IS1111.fasta    ISSm3_IS110_unknown.fast

indir=Input 
outdir=${indir}/Refined_range_fastas
mkdir -p ${outdir}

list=Input/IS110_IDs_with_groups.txt

for fasta in ${indir}/IS*fasta
do
	if [[ "$fasta" =~ .*IS1111.fasta ]]
	then
		# Group is IS1111 sub-group of family IS110, use character range 1000 bp to end of sequence:
		pattern="1000-"
		
	elif [[ "$fasta" =~ .*unknown.fasta ]]		
	then		
		# Group is IS110 family, no sub-group, use character range 1 - 500 bp:
		pattern="1-500"			
	fi
		
	header=$(awk 'NR==1' ${fasta})
	seq=$(sed -e "1d" ${fasta} | tr -d '\n')		
	shorter_seq=$(echo $seq | cut -c $pattern)		
	out=${outdir}/$(basename $fasta | sed "s/.fasta/_${pattern}bp.fasta/")
		
	# Print new fasta
	echo Writing ${out}
	printf "${header}\n${shorter_seq}\n" > ${out}	

done
		
