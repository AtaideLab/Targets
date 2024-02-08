# Revised BLAST filtering parameters

Revisions commenced 23/9/23 by Cali Willet [PIPE-4504](https://ctdshub.atlassian.net/browse/PIPE-4504) 

## Client request

1.	Re-run BlastN with % Identity 95-100%, E value 0, keep rest same.
2.	Re-run BlastN with % Identity 95-100%, E value 10, keep rest same.
3.	Run your extract flanks script with +-60 bp flanks, concatenate.

Name the output folders as: 

- Flanking_Fastas_60bp_IdentI95_E0
- Flanking_Fastas_60bp_Ident95_E10

Once you have the concatenated flanks, that's all I need and I've written some scripts to take care of the rest from here. 

## Compute

NCI Gadi HPC

```
/scratch/er01/PIPE-4504-IS-BLAST
cd /scratch/er01/PIPE-4504-IS-BLAST
git clone https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST.git
cd PIPE3657-IS-BLAST
mv Output Output-2022
mkdir Output
Output-2022/IS110_complete.bacterial.blast.out Output
```

Rename of original output files to avoid any unintentional over-writes. 

Full unfiltered BLAST output was saved so there is no need to re-run BLAST. Copied it from `Output-2022` to `Output` for compatability with workflow.  

## Adjust parameters in BLAST filter script

Script name: `Scripts/filter_first_round_blast.pl`

Original script did not have the filter params coded as variables. Changed these to variables and added a 'filterName' variable to create unique output filenames using the E value and percent identity values:

```

#########################################################
# FILTERING PARAMETERS: TO BE ADJUSTED BY USER
#########################################################

# Subject length rnage. To turn off filtering by subject length, 
# change min_subject_length to zero and change max_subject_length
# to a very large number 
my $min_subject_length = 0; 
my $max_subject_length = 100000000000; 

# Filter for hits with e_value less than or equal to:
my $max_e_value = 0; 

# Filter for hits with percent identity of greater than or equal to:
my $min_pc_ident = 95; 


# Name of filter, to be used to name the output files
# Subject length filter not included in auto-generated filter name, 
# over-ride the filterName variable manually if required
my $filterName = "Ident$min_pc_ident\_E$max_e_value"; 


######################################################### 
```

## Number of filtered BLAST hits

- Email exchange with Rezwan resulted in removal of the subject length filter
- All unfiltered hits have less than 10 E value
- There were 77,324 hits with E value equal to zero
- Once the other filter parameters were applied, the E value zero filter did not make a difference to the output compared to max E value 10
- Confirmed with Rezwan if any other param changes were required, he has asked to move forward with only Ident95_E0
- This gives ~ 8X more hits compared to the original parameters of E value 0, % identity 100, subject length range 100 Kb - 10 Mb 

```
[cew562@gadi-login-05 PIPE3657-IS-BLAST]$ wc -l Output/*filtered
   32317 Output/IS110_complete_Ident95_E0.bacterial.blast.filtered
   32317 Output/IS110_complete_Ident95_E10.bacterial.blast.filtered
    4097 Output/IS110_complete_ORIGINAL.bacterial.blast.filtered
``` 

Of 331 IS with BLAST hits, 259 passed the original filter and 308 have passed the new filters. 
```
[cew562@gadi-login-09 PIPE3657-IS-BLAST]$ awk 'NR>1  && $5>0' Output/IS110_complete_Ident95_E0.bacterial.blast.report | wc -l
308
```


## Rerun flank extraction

Script name: `Scripts/extract_flanks_keepNonUnique.pl`

This script extracts 200 bp flanks from the subject sequences from the BLAST-formatted nr/nt database. The flank extraction uses BLAST+ `blastdbcmd` utility. It takes around 6 seconds per entry (6.5 hours for the original 4097 filtered hits).
Given there are now 32,317 hits, it is worthwhile to make this script parallel.

Writing new script `Scripts/extract_flanks_keepNonUnique_parallel.pl` 

Reconfiguring the workflow to be parallel. This parallelisation won't require changes to the main part of the script, which writes a multi-fasta per IS, just a change to the input file names and a pre-step to make the parallel inputs.

Writing new scripts:
```
Scripts/extract_flanks_parallel_make_inputs.pl
Scripts/extract_flanks_keepNonUnique_parallel.pl
Scripts/extract_flanks_keepNonUnique_parallel.pbs
```

### Make parallel inputs file 

`Scripts/extract_flanks_parallel_make_inputs.pl` is run on the login node. It reads the filtered BLAST output and makes N new output files, holding all of the hits but in shorter lists. N is equal to the number of chunks specified by the user. 
I have used a chunk value of 224, being suited for 8 Broadwell nodes. Each chunk for this dataset has 145 hits and should complete in < 20 minutes. Due to rounding, there may be 1 file less than the number of chunk files specified. This will not harm overall efficiency very much for large chunk numbers.  

Open the script and edit variables for `chunks` (number of parallel chunks of tasks to run), `filterName` (name of the filter emit by the previous BLAST filter step) and `dataset` (name of the dataset, used throughout the workflow). 

To run on the login node:
```
bash Scripts/extract_flanks_parallel_make_inputs.pl
```

Output fill be BLAST filter outputs but spread of <chunks> number of files, and a TEMPLIST that is used to parallelise the tasks. 

### Run parallel flank extraction 

Parallelisation is by `nci-parallel` utility within the script `Scripts/extract_flanks_keepNonUnique_parallel.pbs`, which allocates the N chunks of hits across multiple nodes and runs `Scripts/extract_flanks_keepNonUnique_parallel.pl` in paralel, once per shorter input hit list.

To run, update resources for he number of chunks and expected walltime, allowing 6 seconds per hit. For example, for 223 chunk files (requested 224 chunks to script but rounding led to one less file) and 32,216 hits input, this is 145 hits per file so 145 X 6 seconds ~ 14.5 minutes. Allow 1 Broadwell CPU per chunk/task and 9 GB RAM per chunk/task. Use all the JOBFS on the node, so 8 nodes X 400 GB jobfs = 3200 GB jobfs. So resources for this job would be 224 CPU, 2016 GB RAM, 3200 GB jobfs, and allow 1 hour waltime "just in case". Jobfs is required as there are a large number of small files made at once and this affects performance of the Lustre filesystem.

The script assings a "-<chunkNum>" suffix to the multifastas, and these are then merged per IS at the next step. There will be no unintentional filename clashes during parallel run because the temp fastas are named with a combination of IS ID, subject accession, subject start and subject end. To check this:

```
awk -F "\t" 'NR>1 {print $1"_"$8"_"$10"_"$11}' Output/IS110_complete_Ident95_E0.bacterial.blast.filtered | sort | uniq | wc -l
#32316
```


Save then sumbit with:
```
qsub Scripts/extract_flanks_keepNonUnique_parallel.pbs
``` 
The number of output fastas should be equal to the number of lines in the TEMPLIST inputs file. The total line count of these fastas should equal the total number of hits (in this case 32,316) times 2 (for fasta headers).  

#### Unusual run times

When testing a single chunk of 145 files, the run time was reliably 14-15 minutes, ie 6 seconds per hit. 

When the whole 223 chunks were run together over 8 nodes, the run time per hit was more like 20 seconds. 

Have given 1 hour to the job but I expect it to fail...has done 3567 hits of 32 K in  26 minutes... This extrapolates to 236 minutes. I wonder if its run time is affected by multiple reads on the database? In which case, parallelising to a greater number of chunks would not help. I now recall some strnage run time behavioru when attempting to benchmakr the initial BLAST. 

Killed it and resubmitted with 12 hrs walltime. It completed in 3 hrs 51. 

#### Check the output:
```
$ grep "exited with status 0" PBS_logs/extract_flanks_parallel.e  | wc -l
223
$ wc -l ./Output/IS110_complete_Ident95_E0.bacterial.blast-TEMPLIST
223 ./Output/IS110_complete_Ident95_E0.bacterial.blast-TEMPLIST
$ cd  Output/Flanking_fastas_Ident95_E0/
$ grep DONE *log | wc -l
223
$ wc -l *fasta-*
 64632 total
```

Number of fasta lines equals expected (32,316 X 2) and no errors detected. 


### Combine parallel outputs per IS

Take all the 400 bp flank fastas per IS and concatenate them into one multi-fasta per IS sequence. 

Adjust `filterName` and `dataset` variables, then run:

```
bash Scripts/combine_flanks_from_parallel.sh
```

#### Check the output: 

The number of output fastas in `./Output/Flanking_fastas_<filterName>` should equal the number of IS with >0 hits passing the previous BLAST filter:
```
$ cd Output/Flanking_fastas_Ident95_E0
$ awk 'NR>1 && $5>0' ../IS110_complete_Ident95_E0.bacterial.blast.report | wc -l
308
$ ls -1 *fasta | wc -l
308
```

The total line count of these fastas should equal the number of inputs (in this case, 32,316) times 2 (for fasta header):

```
$ wc -l *fasta
64632 total
```

Once this step has been completed and checked, the temp chunk fastas and status logs can be deleted:
```
\rm -rf ./Output/Flanking_fastas_<filterName>/*_400bp_flanks.fasta-*
\rm -rf ./Output/Flanking_fastas_<filterName>/*log
```

 
### Clean up temp files

After the parallel outputs have been created and combined, all TEMP files in the output dir can be deleted:

``` 
\rm -rf ./Output/*TEMP*
```
 
## Extract 60 bp flanks

From the 400 bp flank output, extract a shorter region of 60 bp flank. 

Script: `Scripts/extract_shorter_flanks.pl`

This was returning many 'substr out of range' errors, because we have now dropped the minimum subject length filter. Sequences with very few or even zero flanking base pairs have been let through, for example:
```
#qseqid qlen    length  qstart  qend    sseqid  stitle  sacc    slen    sstart  send    pident  mismatch        gapopen evalue  bitscore
IS1111C_IS110_IS1111    1452    1452    1       1452    gi|145004|gb|M80806.1|COXTRANSPO        Coxiella burnetii transposase (IS1111a) gene, complete cds      M80806  1450    1       1450    99.036  12      2       0.0  2603
```

Edited the script to allow for left and right flanks that are fewer than 200 bp, however I now have concern that a small number of the 60 bp flanks previously generated may not be spot on in terms of the positions given to the perl 'substr' function. For IS that were within 200 bp of the end or start of the subject sequence, the 60 bp flank will be wrong, because they were all extracted from 140 bp position onwards. The extract flanks script that used `blastdbcmd` tool allowed for this, but I failed to add the allowance into the `extract_shorter_flanks` scripts. I expect the number of sequences affected will be small, but I will rerun anyway and advise Rezwan that this has been done and which fasta have been affected. 


Emailed Rezwan to ask if he would like to apply a size filter to the flanks: yes he would only like to keep full-size flanks, so 120 bp in this case. Also asked if he would like the final steps of MSA and Emboss Cons from the previous workflow to be run; he does not require this. 

The script takes input flank size as first and only command line argument, and filters out all input sequences that are < 2 X flank size and filters out all output flank sequences that are less than 2 X flank size. Failed sequence headers are recorded. 

To run (a few seconds on the login node):
```
$ time perl Scripts/extract_shorter_flanks.pl 60
Total input sequences: 32316
Total 120 bp output flank sequences: 32244
Total failing input length filter of 120 bp: 41
Total failing output length filter of 120 bp: 31

Failed sequence headers are written to file ./Output/Flanking_fastas_60bp_Ident95_E0/IS110_complete_60bp_Ident95_E0.failed.txt

New 120 bp flank fastas are written to directory ./Output/Flanking_fastas_60bp_Ident95_E0
```

### Check for flank length on original output

Copied the above script and edited to check if any o the original 60 bp flanks output was affected by the IS insertion site residing within 200 bp of the subject start or end. 

Only one sequence was affected. Emailed Rezwan to advise.
```
$ perl check_extract_original_shorter_flanks.pl 60
Total input sequences: 4096
Total 120 bp output flank sequences: 4095
Total failing input length filter of 120 bp: 0
Total failing output length filter of 120 bp: 1

Failed sequence headers are written to file ./Output-2022/Flanking_fastas_60bp_corrected/IS110_complete_60bp.failed.txt

New 120 bp flank fastas are written to directory ./Output-2022/Flanking_fastas_60bp_corrected

$ cat ./Output-2022/Flanking_fastas_60bp_corrected/IS110_complete_60bp.failed.txt
WARN OUTPUT: New flank sequence for >ISEc32_IS110_unknown_within_CP088663_leftFlank_1-15_rightFlank_1196-1395 less than target length 120: length 75 
```
 
Deleted the affected sequence from within the flank fasta for that IS at his request:
  
```
# Line count shows the corrected file has one fewer sequence, as expected:
$ wc -l Flanking_fastas_60bp/ISEc32_IS110_unknown_60bp_flanks.fasta 
264 Flanking_fastas_60bp/ISEc32_IS110_unknown_60bp_flanks.fasta
$ wc -l Flanking_fastas_60bp_corrected/ISEc32_IS110_unknown_60bp_flanks.fasta 
262 Flanking_fastas_60bp_corrected/ISEc32_IS110_unknown_60bp_flanks.fasta

# Check one random other IS to confirm no difference in the rest of the output:
$ sdiff -s Flanking_fastas_60bp/ISYen1_IS110_IS1111_60bp_flanks.fasta Flanking_fastas_60bp_corrected/ISYen1_IS110_IS1111_60bp_flanks.fasta

# Replace the affected fasta with the corrected one:
$ mv Flanking_fastas_60bp_corrected/ISEc32_IS110_unknown_60bp_flanks.fasta Flanking_fastas_60bp

# Delete the remaining unaffected fastas:
$ \rm -rf Flanking_fastas_60bp_corrected/
```










































































