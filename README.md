# Workflow to generate conserved target sites of IS1111 and IS110 family of insertion sequences

## Overview 

Insertion sequences (IS) transpose into genomes in a sequence specific manner. This repository contains a collection of scripts to find the target sites where the IS transposes into. This workflow is implemented in the Linux/Unix command-line environment.

A multi-fasta containing all query IS sequences is first searched against the bacterial and archaeal sequences within the NCB non-redundant nucleotide database using `blastn` from the `BLAST+` package. The BLAST output is then filtered for E value 0 and minimum sequence identity 95%. The `blastdbcmd` tool is then used to extract 200 bp of flanking seqeunce from either side of the insertion sequence. All sequences including exact duplicates are retained in a multi-fasta of insertion flanks per IS. From the 400 bp concatenated flank fastas, 20 bp flanks are extracted, taking care to ensure that the junction of the left and right flank is preserved in the centre of the 40 bp output sequence. These fastas are then processed with `WebLogo`. The resultant image files can then be viewed to identify insertion sequence recognition sequences.

The workflow here focuses on two IS families: IS110 and IS111. We describe the creation of a multi-fasta for these families from the [ISFinder](https://isfinder.biotoul.fr/) database and include relevant inputs and outputs [REZWAN/SANDRO - DO YOU WANT TO INCLUDE OUTPUTS?]. In order to replicate this workflow on your chosen IS family, a multi-fasta of IS sequences is required. Ensure to update relevant filepaths in the scripts to match our input dataset.  

A demo/test dataset is provided with two IS, one from each of IS110 and IS1111. Results can be compared with ours published here to check the workflow before implementing on your dataset. 

## Compute requirements
The demo workflow has been tested on a Linux command-line environment. The full workflow was implemented on [NCI Gadi HPC](https://nci.org.au/our-systems/hpc-systems) running PBS Pro on CentOS. 

## Software dependencies
- BLAST+
- perl
- python
- python packages `biopython`, `bio`, `weblogo`

## Running the demo
1) Run the BLAST against mini db
Notes: this script requires `blast+` module, it includes a `module load blast+` command.  If `blast+` is already in your path, you can delete/hash out this line, or edit to suit the requirements of your environment.  

Change into the base working directory `demo`, then run:

```
$ bash Scripts/blast_IS.pbs test
Running two IS test
```

Output:
```
$ wc -l Output/IS_2sequence_demo.bacterial_archaeal.blast.out 
123 Output/IS_2sequence_demo.bacterial_archaeal.blast.out
```

2) Filter the BLAST for minimum identity 95% and E value 0:

```
$ perl Scripts/filter_blast.pl
```

Output:
```
$ wc -l Output/IS_2sequence_demo_Ident95_E0.bacterial_archaeal.blast.filtered 
122 Output/IS_2sequence_demo_Ident95_E0.bacterial_archaeal.blast.filtered
$ cat Output/IS_2sequence_demo_Ident95_E0.bacterial_archaeal.blast.report 
#Sequence_ID    Family  Group   Raw_hits        Passing_hits
ISPlge4 IS110   IS1111  67      67
ISPsy35 IS110   NA      56      54
```
2 hits failed filtering. 

3) Create flank span/range lists for batch flank extraction:

```
$ perl Scripts/extract_flank_ranges.pl
```

Output:
```
$ wc -l Output/Flanking_fastas_Ident95_E0/*
    0 Output/Flanking_fastas_Ident95_E0/failing_flank_warnings.txt
  121 Output/Flanking_fastas_Ident95_E0/left_flank_ranges.batch.txt
  121 Output/Flanking_fastas_Ident95_E0/right_flank_ranges.batch.txt
```

4) Extract 200 bp flanking sequence for hits in filtered BLAST output using `blastdbcmd`:

```
$ bash Scripts/extract_flanks_submit.sh test
Running 2 IS demo: creating left flanks
Running 2 IS demo: creating right flanks
```

Output:
```
$ wc -l Output/Flanking_fastas_Ident95_E0/*fasta
  242 Output/Flanking_fastas_Ident95_E0/left_flanks.fasta
  242 Output/Flanking_fastas_Ident95_E0/right_flanks.fasta
```


5) Concatenate the flanks into one multi-fasta per IS:
```
$ perl Scripts/extracted_flanks_postprocess.pl
```

Output:
```
$ wc -l Output/Flanking_fastas_Ident95_E0/200bp_flanks/*
  134 Output/Flanking_fastas_Ident95_E0/200bp_flanks/ISPlge4_IS110_IS1111_200bp_flanks.fasta
  108 Output/Flanking_fastas_Ident95_E0/200bp_flanks/ISPsy35_IS110_unknown_200bp_flanks.fasta
```

6) Extract smaller flanks of desired length for WebLogo generation

Provide desired flank length (eg 20 bp) as first and only command line argument. Flanks failing the size filter will be printed to a list file in the output directory. 

```
$ perl Scripts/extract_shorter_flanks.pl 20
Total input sequences: 121
Total 2 x 20 bp output flank sequences: 121
Total failing input length filter of 40 bp: 0
Total failing output length filter of 40 bp: 0

Failed sequence headers are written to file ./Output/Flanking_fastas_Ident95_E0/20bp_flanks/target_length_failed.txt

New 40 bp fastas are written to directory ./Output/Flanking_fastas_Ident95_E0/20bp_flanks
```

7) Create WebLogos

This step requires `biopython`, `bio` and `weblogo` python packages. If you do not have these installed, run:

```
pip install biopython
pip install bio
pip install weblogo
```

Provide the directory containing the fasta you wish to create WebLogos for as a command-line argument to the script. Note that this script requires all input sequences per IS are of equal length. The downstream script `extract_shorter_flanks.pl` ensures this. The 'filter_name' and 'flank_size' component of the input directory path are used to define the output directory path.

```
$ python3 Scripts/weblogo_multipng.py Output/Flanking_fastas_Ident95_E0/20bp_flanks/
Creating WebLogos on fastas in Output/Flanking_fastas_Ident95_E0/20bp_flanks/
Writing WebLogos to Output/WebLogos/Ident95_E0_20bp_flanks
Processing: ISPsy35_IS110_unknown_20bp_flanks.fasta
Processing: ISPlge4_IS110_IS1111_20bp_flanks.fasta
```

8) Compare the output
Your generated output will be in `demo/Output`. The expected output is in `demo/expected_output`. 


## Running the workflow

### Preparing the multi-fasta
- Steps to get all the IS110 and IS1111
    - one script frm the multifasta craetion, the rest is bash one liners: update_missing_sequences_fasta_headers.pl

### BLAST and filter    
blast_IS.pbs 
filter_blast.pl

details TBA

### Extract insertion site flanking sequence
extract_flank_ranges.pl
extract_flanks_submit.sh
    - submitted by above: extract_flanks.pbs
extracted_flanks_postprocess.pbs
    - submitted by above: extracted_flanks_postprocess.pl
extract_shorter_flanks.pl <N> 

details TBA

### Create WebLogos
weblogo_multipng.py



## Citations
- Siguier P et al. (2006) [ISfinder: the reference centre for bacterial insertion sequences](https://pubmed.ncbi.nlm.nih.gov/16381877/) Nucleic Acids Res. 34:D32-D36
- [ISFinder database](http://www-is.biotoul.fr)
- BLAST non-redundant nucleotide database: Sayers E et al. [Database resources of the National Center for Biotechnology Information](https://pubmed.ncbi.nlm.nih.gov/33095870/) Nucleic Acids Res. 49(D1):D10-D17
- Camacho C et al. (2008) [BLAST+: architecture and applications](https://pubmed.ncbi.nlm.nih.gov/20003500/) BMC Bioinformatics 10:421
- Crooks G et al. (2004) [WebLogo: a sequence logo generator](https://pubmed.ncbi.nlm.nih.gov/15173120/) Genome Res. 14(6):1188-90



# Acknowledgements
The authors acknowledge the technical assistance provided by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney and the Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia. The authors acknowledge the use of the National Computational Infrastructure (NCI) supported by the Australian Government and the Sydney Informatics Hub HPC Allocation Scheme, supported by the Deputy Vice-Chancellor (Research), University of Sydney and the ARC LIEF, 2019: Smith, Muller, Thornber et al., Sustaining and strengthening merit-based access to National Computational Infrastructure (LE190100021).


