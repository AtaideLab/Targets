# < rezwan insert your catchy workflow title here>

## Overview 

Insertion sequences (IS) transpose into genomes in a sequence specific manner. This repository contains a collection of scripts to find the target sites where the IS transposes into. This workflow is implemented in the Linux/Unix command-line environment.

A multi-fasta containing all query IS sequences is first searched against the bacterial and archaeal sequences within the NCB non-redundant nucleotide database using `blastn` from the `BLAST+` package. The BLAST output is then filtered for E value 0 and minimum sequence identity 95%. The `blastdbcmd` tool is then used to extract 200 bp of flanking seqeunce from either side of the insertion sequence. All sequences including exact duplicates are retained in a multi-fasta of insertion flanks per IS. These fastas are then processed with `WebLogo`, restricting the logo to the central 40 bp surrounding the insertion point. The resultant image files can then be viewed to identify insertion sequence recognition sequences.

The workflow here focuses on two IS families: IS110 and IS111. We describe the creation of a multi-fasta for these families from the [ISFinder](https://isfinder.biotoul.fr/) database and include relevant inputs and outputs. In order to replicate this workflow on your chosen IS family, a multi-fasta of IS sequences is required. Ensure to update relevant filepaths in the scripts to match our input dataset.  

A test dataset is provided with a small number of IS. Results can be compared with ours published here to check the workflow before implementing on your dataset. 



## Software dependencies
- BLAST+
- perl
- python



## Compute requirements
The IS110/IS1111 family analysis was computed on [NCI Gadi HPC](https://nci.org.au/our-systems/hpc-systems) running PBS Pro. The extraction of insertion sequence flanks was computationally intensive for this dataset (35,520 filtered BLAST hits) and was parallelised to 119 CPU to achieve a walltime of <N>. 



## Running the test

Clone the repository and change into the test directory:
```
git clone https://github.com/AtaideLab/Targets.git
cd Targets/test
```

1) run the BLAST against mini db
Note that this script requires blast+ module, it incldues a 'module load blast+' command. If blast+ is already in your path, can delete/hash out this line. The full versin of this script is a PBS job, uses mt-mode 1, and searches bacteria and archaea from the full non-redundant nucleotide database. This demo script has omitted those features for simplicity of executing the test run. The test database has been restricted to include only the accessions that the 2 test IS align to.  
Ignore the PBs component, just run as bash: 
```
$bash bin/blast_IS.pbs
```

Creates output: `Output/IS_2sequence_demo.bacterial_archaeal.blast.out`

2) filter the BLAST
```
$ perl bin/filter_blast.pl
```

Creates output: `Output/IS_2sequence_demo_Ident95_E0.bacterial_archaeal.blast.report` and `Output/IS_2sequence_demo_Ident95_E0.bacterial_archaeal.blast.filtered`

3) create ranges files for batch flank extraction

```
$ perl bin/extract_flank_ranges.pl

$ ls -1 Output/Flanking_fastas_Ident95_E0/
right_flank_ranges.batch.txt
left_flank_ranges.batch.txt
failing_flank_warnings.txt
```

4) extract 200 bp flanks

This script usually submits 2 jobs to PBS, but this part has been hashed out for the test run to just submit the 2 jobs in series on the login node. It still completes in < 10 seconds. 
```
bash bin/extract_flanks_submit.sh
```

Creates 2 new files in  `Output/Flanking_fastas_Ident95_E0` directory: `left_flanks.fasta` and 
`right_flanks.fasta`. 
 

5) concatenate flanks
For full workflow, this is launched by PBS wrapper `extracted_flanks_postprocess.pbs`. Not required here - 10 seconds to run on login node. 

```
perl bin/extracted_flanks_postprocess.pl
```

New outputs are 400 bp fastas in `Output/Flanking_fastas_Ident95_E0/200bp_flanks`, one multi-fasta per IS: 

```
$ ls -1 Output/Flanking_fastas_Ident95_E0/200bp_flanks/
ISKpn60_IS110_unknown_200bp_flanks.fasta
ISPa11_IS110_IS1111_200bp_flanks.fasta
```

6) extract smaller flanks of desired length for WebLogo generation

Provide desired flank length (eg 20 bp ) as first and only command line argument:

```
$ perl bin/extract_shorter_flanks.pl 20
Total input sequences: 1691
Total 2 x 20 bp output flank sequences: 1690
Total failing input length filter of 40 bp: 0
Total failing output length filter of 40 bp: 1

Failed sequence headers are written to file ./Output/Flanking_fastas_Ident95_E0/20bp_flanks/target_length_failed.txt

New 40 bp fastas are written to directory ./Output/Flanking_fastas_Ident95_E0/20bp_flanks
```


7) Create WebLogos

TBA



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


