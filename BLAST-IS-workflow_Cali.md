---
title: "PIPE-3657: IS BLAST search query"
author: "Cali Willet"
date: "2022-11-25"
output: html_document
---

# Research group
## Lead chief investigator
Dr Sandro Fernandes Ataide <sandro.ataide@sydney.edu.au>
Senior Lecturer in Structural Biology, SoLES

## Research Associate
Rezwan Siddiquee <rezwan.siddiquee@sydney.edu.au>

# Project overview

Find the conserved target sites of IS110 and IS1111

# Data 

Data is all publicly available. 

The IS seqeunces are available on the [ISFinder website](https://isfinder.biotoul.fr/). This website does not allow bulk data downloads, although there is a [github repository](https://github.com/thanhleviet/ISfinder-sequences) from 2020 with a fasta data dump.
Dr Ataide has contacted ISFinder group for bulk download permission, as per the request on the ISFinder website. 

The 350 IS sequences from the two families will be searched against the non-redundant nucleotide database from NCBI. 

# Analysis overview

* Input: 350 IS from [ISFinder](https://isfinder.biotoul.fr/) 
    * IS110 family
    * IS1111 subgroup of family IS110  
* Develop script(s) that performs for each IS:
    * Nucleotide BLAST against the standard nrnt database
    * Filter matches to some cut-off criteria (TBD)
    * Obtains 200 bp flanking sequence, concatenates the flanking sequence into 400 bp string

# Analysis considerations and challenges

For the first blast of IS against nr, this must be stringent - eg > 95% identity for full length of IS. 

Some genomes have multiple copies of the same IS, and thee may be found int he 5'->3' or 3'_.5' orientation. For the ClsutalW alignment, some sequecnes will need to be flipped. Is there a way to determine this part via code, or will it need to be assessed by eye and then re-run that modular part of code to do the flipping and re-alignment?

Rezwan has sent an example excel sheet showing the manual process of finding hits and flanks. Use these 2 examples(one from each group) as starting sequences to re-create the process computationally. 

Use BEDtools get fasta to extract the flanking sequence from the nr database, eg `bedtools getfasta -fi sequence.fasta -bed test.bed`

# Analysis commencement 25/11/22

# Obtain the 350 sequences in fasta

Obtain seqs from github db, although the ID names are not matching what they are supposed to based on the readme - compare the list of IDs to those 350 from ISFinder provided by Rezwan. Rezwan to get the missing (newer) sequences if not too many.

 `awk '$1~/^>/' ./ISfinder-sequences/IS.fna | sed 's/^>//' | grep _IS110_ | wc -l`

There are 331, so only 19 sequences have been added since this data dump in 2020. 
Edit 2/12/22: 331 was taken from the repo I had used for a project in 2020. After I cloned the repo again, there are now 335 sequences. These 4 (ISEfa16 ISKpn60 ISKqu3 ISLad6 ) are already covered in the 19 squences Rezwan has provided in multi-fasta, so I will use the newer repo fasta and remove these 4 manually from the file Rezwan sent

Clone the repo with IS fasta:
`git clone https://github.com/thanhleviet/ISfinder-sequences.git`

Make a list of IDs: 

`awk '$1~/^>/' ./ISfinder-sequences/IS.fna | sed 's/^>//' | grep _IS110_ | cut -d '_' -f 1 > IS110_IDs_from_github.txt`

Find the sequences missing from Github repo:

`comm -13 IS110_IDs_from_github.txt-sorted IS110_family_ID_list.txt-sorted > IS110_19seqs_missing_IDs.txt`

There are 19 [ see edit note above, updated number is 15 ], emailed these to Rezwan and requested fasta or multi-fasta that I can combine with the other 331 [335]. 

Extract the 331 [335] IS110 members from the Github repo to a multifasta: 
 ```
awk '$1~/^>/' ./ISfinder-sequences/IS.fna | grep _IS110_ | sed 's/^>//' > IS110_fastaHeaders_from_github.txt
module load seqtk/1.3
seqtk subseq ./ISfinder-sequences/IS.fna IS110_fastaHeaders_from_github.txt > IS110_from_Github.fasta
```

Concatenate the 19 [ 14] manually obtained by Rezwan with the 331 [ 335 ] from Github:
Rezwan provided multi-fata for 18 sequences, as one sequence (IS1111B) was missing from ISFinder.
Convert to unix line endings, remove blank lines:
```
dos2unix Multi_Fasta_of_Missing_IS_from_IS110_IS1111.fasta
sed -i '/^$/d' Multi_Fasta_of_Missing_IS_from_IS110_IS1111.fasta
```

Open the fasta, manually remove the 4 duplicated sequences, then add the family ID to the IS name, to be consistent with the fasta headers from the Github repo. The Github fasta has the following fasta header format:

<ID>_<family>_<group>

If there is no group, then 'unknown' is used. For the 14 sequences remaining in the fasta provided by Rezwan, add the additional header info  with perl script (should have asked Rezwan to do this initially). Check, then replace the original fasta, then concatenate the 14-seq fasta with the 335-seq fasta:

```
perl update_missing_sequences_fasta_headers.pl
mv Multi_Fasta_of_Missing_IS_from_IS110_IS1111.fasta-reheadered Multi_Fasta_of_Missing_IS_from_IS110_IS1111.fasta
cat IS110_from_Github.fasta Multi_Fasta_of_Missing_IS_from_IS110_IS1111.fasta > IS110_complete.fasta
```

There is now one fasta with all 349 available IS110 family sequences: `IS110_complete.fasta`

# Compute
Goal is to use Gadi to perform the analysis, as Gadi is faster, however the scripts will be portable to Artemis (apart from a few small directive changes) 

# BLAST IS against NCBI nr/nt database

BLAST on gadi:
`module load blast+/2.13.0`

02-12-22: updating the NCBi preformatted blast databases, to `/g/data/er01/NCBI/preformatted_2022_11_21`

Working directory is `/scratch/er01/cew562/PIPE-3657-IS`


## Blast one sequence

IS1111 member ISPa11, example seqeunce from Rezwan. Obtain one blast output to setup the scripts with.

Target seq: GTAGGGCGAATAAC

Rezwan: "I want to verify this seq and also find an extended version.Manually BlastN and find a few copies"



If you need FASTA for selected sequence(s) from these BLAST databases, you can obtain it as follows (the sequence of interest is identified by the accession u00001 in this example):

`blastdbcmd -entry u00001 -db nr -out u00001.fsa`


Using blastn output format 6, which has the following headers (not emitted with blastn output):

`qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`

To create custom output format, see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6 

First round BLAST params from Sandro and Rezwan:
* E-value of 0
* %identity 100% (I've found recently that any less than 100% is not useful for us as the target won't be accurate enough for our downstream experiments). 
* Alignment length - not sure what you mean by this. Won't the alignment length vary because each "IS" is a different length+200 bp flanks?
* 100% query coverage
* Another thing to note is to check the accession length to make sure it looks like a bacterial genome size like 4 x10^6 bp. Sometimes you'll find 100% match with a much smaller accession length because those are independent integrons/plasmids

By applying the stricter filters aBbove, we get 97 perfect hits for that sequence. Would you be wanting to obtain and BLAST the flanking sequence for all 97 of those? Is it likely that some of the 97 flanking sequence fastas are identical ( in which case we would not need to blast these, that’s what I meant by unique flank fastas)? 

Reply from Rezwan: "We would want to obtain those 97 flanks, then run a Clustal (or any other) alignment on them to find out the boundaries upto which they are conserved. 
So we start with 200+200 = 400, then only 50 bp on either side might be conserved. 
Then we combine that 50 bp and blastn again, to find out if some bases have been inserted in between or is it perfect. 

Using the following headers for custom blast output format:
`qseqid qlen length qstart qend sseqid stitle sacc slen sstart send  pident mismatch gapopen evalue bitscore`


## Small amount of benchmarking on Gadi
1. Run ISPa11 blast on one whole normal node
2. Run ISPa11 blast on 1/2 normal node
3. Run ISPa11 blast on 1/2 hugemem node
4. Run ISPa11 blast on 1/4 hugemem node
5. Repeat the above steps for a 6-sequence fasta

The % CPU and mem is VERY LOW!

Identical run time for 12 or 24 CPU on hugemem. 
Faster for 6 seqs on hugemem than for 1. 
Could there have been a pre-load of index into RAM, and the 1seq and 6seq jobs shared the same node? Yes this appears to be possible, see resource excel sheet for shared nodes on 2 jobs (inc the fastest one) 

## Make the search refined to bacteria
https://www.biostars.org/p/420614/

2 is the tax id for order bacteria 

Downloaded edirect to gdata apps and added path to blast+ .base file
`sh get_species_taxids.sh -t 2 > bacterial.taxids`

Added `taxidlist` flag to blast command, resubmitted 2 jobs - the slowest (1 seq 48 CPU normal) and the fastest (6 seqs 12 CPU hugemem) 

Also to explore:
1. Different threading option ["threading by query"](https://www.ncbi.nlm.nih.gov/books/NBK571452/) - although the specs of this job do not really fit with the recommendation of bases per thread
2. mpirun - according to [here](https://hpc.ncsu.edu/Software/Apps.php?app=BLAST) just adding the mpirun part, but given there are specific mpi-built blast tools, i doubt its that simple... 
3. nci-parallel, running 48 x 1 CPU tasks on hugemem (but this is not portable to Artemis)

Tested 1, there was a warning in the .e log file yet the results were still faster. Expect ~ 14 hrs for the 349 seqs on one hugemem node

## Run the first round BLAST of all IS against the 'bacterial' taxids in nt/nr

- This is "Step 1" from [project steps](https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST/blob/master/BLAST-IS-workflow.md#project-steps-laid-out)

-Run on Gadi:
`qsub Scripts/blast_IS_complete.pbs` 

Wow this job completed in 12 minutes! How bizarre, blast is so unpredictable! 

## Filter the IS110_complete blast output 

- Also part of  "Step 1" from [project steps](https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST/blob/master/BLAST-IS-workflow.md#project-steps-laid-out)

Filter the initial BLAST on the above supplied parameters, and create a report of the BLAST results

Parameters applied: 
•	E value of 0
•	% identity 100.000 
•	Subject length >= 4,000,000 (also restricted the BLAST search to bacterial taxonomic IDs within the non-redundant nucleotide database, given that group is only interested in bacterial genome hits)
•	Query coverage 100%

Run: 
`perl Scripts/filter_first_round_blast.pl`

Creates: 
* IS110_complete.bacterial.blast.filtered
    * BLAST hits passing filters, with the addition of a header line 

* IS110_complete.bacterial.blast.report
    * Sequence ID, group status, raw and filtered hit counts
    * 331 of the 349 sequences have results, range of hits is from 1 to 39,124
    * 211 of the sequences have zero passing hits
    * 5/12/22: emailed this report to Sandro and Rezwan to find out if they want relaxed filters for these zero-hit sequences
        * Changed filters to min 100 Kbp subject size, and max 10 Mbp
            *  Now only 90 seqs with no passing hits
        * Manual example of Kpn4 provided
            * Question regarding whether or not the 3 to 5 orientation is present, as Kpn4 should have more than 1 copy on 'Klebsiella pneumoniae strain RGF99-1 plasmid', which Rezwan notes are mostly in 3 to 5 
            * Blast checks for both orientations by default (-strand parameter, default = both) 
            * Web blast shows only 1 Kpn4 hit to RGF99-1, same as cmdline
            * Web shows 'plus/minus' hit and cmdline describes this by having qstart > qend
            * 6/12/22 Sandro confirmed Kpn4 has only 1 insertion; noted EC11 and Pa11 as IS with multiple insertions. Use these later  to check the scripts are finding multiple insertions properly. 

From https://doctorlib.info/medical/blast/7.html 
"The plus strand is the sequence in the FASTA file. The minus strand is the reverse complement of this sequence. If the similarity between the query and subject sequences is on the same strand, both sequences are labeled as being on the plus strand and the coordinates increase from left to right. When the minus strand of the query sequence is similar to a database sequence, In NCBI-BLAST, the database sequences are flipped"

## Obtain flanking sequence to create target site multi-fasta

- Last part of  "Step 1" from [project steps](https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST/blob/master/BLAST-IS-workflow.md#project-steps-laid-out)


Need too use NCBI tool to extract fasta from the pre-formatted BLAST database. 
` blastdbcmd -entry CP034908.2 -db $nt -out test.hit -range 50-100`


*  Flank extraction perl script:
`Scripts/extract_flanks.pl` to be executed on compute node with `Scripts/extract_flanks.pbs`
- 16/12/22 
    - Some of the 400 bp flanks are duplicates per IS
     - Emailed questions to Rezwan and Sandro:
         - how to handle this, eg do we need to make the flanks longer?
        - Keep all headers in the output, or just keep the first header for dup/identical flank seqs and record the counts.
            -  Sandro replied to keep the first only and record count

- 12/1/23
    - Script is now working perfectly. 
    - Will need compute node, test took ~ 3 mins 17 for 72 hits 
    - Input file `IS110_complete.bacterial.blast.filtered` has 4097 lines
        - Expect ~ 3-3.5 hrs for full run on one CPU --> quick enough to not warrant splitting the input

 - Run the script on entire blast filtered input, use compute node 

 ```
 qsub Scripts/extract_flanks.pbs
 ```

 - 13/1/23 completed in 6hrs 49 mins. Used all of RAM - I wonder if it will run faster ona hugemem CPU
- There is one error message from blast in the .e log "Error: [blastdbcmd] Failed to parse sequence range (start cannot be empty)"
        - Should have included a check for this:

```
[cew562@gadi-login-07 PIPE3657-IS-BLAST]$ perl Scripts/checks_spans.pl 
WARN: Line num 990 right flank end is 142242 and subject length is 142200
WARN: Line num 1057 left flank start is -184
```
    - Update extract_flanks.pl to correct for this and re-run. 
    - Submittedon hugemem 13/1/23

17/1/23 Script completed with no error now that there is a check to make sure the flanks do not exceed the bounds of the subject sequence 


## Summary of scripts in this pipeline 
- I have not included the preparatory steps that involved getting the input multifasta in order from the 2 online sources.

1. Gadi HPC:
    ```
    qsub Scripts/blast_IS_complete.pbs
    ```
2. Login node:
    ```
    perl Scripts/filter_first_round_blast.pl
    ```
3. Gadi HPC:
    - a. The initial 6.5 hour run t hat used blastdbcmd to extract fasta from BLAST-formtatted db
    ```

    qsub Scripts/extract_flanks.pbs # runs extract_flanks.pl - updated to extract_flanks_keepNonUnique.pl 

    ```
    - b. Script to adjust flank size (up to 200 bp) and ad redundancy for non-unique flanks
    ```
    perl Scripts/extract_shorter_flanks_with_redundancy.pl
    ```





