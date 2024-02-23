#!/bin/bash

set -e 

# Submit compute job for each list of BLAST hits made by
# extract_flanks_make_chunks.pl

# Split lists of BLAST hits to get flanks from
chunk_dir=Output/extract_chunk_tmpdir

# Script to run:
script=Scripts/extract_flanks_parallel.pl

compute_logs=PBS_logs/extract_flanks
mkdir -p $compute_logs

module load blast+
blastdbcmd=$(which blastdbcmd)

project=xh27
storage='scratch/er01+gdata/er01'

for file in ${chunk_dir}/extract_chunk_*
do

	echo $file
	name=$(basename $file | sed 's/extract_//')
	qsub -P ${project} \
		-l storage=${storage} \
		-N ${name} \
		-q normalbw \
		-l ncpus=1 \
		-l mem=4GB \
		-l walltime=02:00:00 \
		-l wd \
		-W umask=022 \
		-o ${compute_logs}/${name}.o \
		-e ${compute_logs}/${name}.e \
		-v filtered_hits="${file}",blastdbcmd="${blastdbcmd}" \
		${script}
	
	sleep 2
done
