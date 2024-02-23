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

-  Add 'tree' of contents
- add steps
- add instruction re comparing outputs (comapre the logos created for the test IS to those in the main outputs dir)



## Running the workflow

### Preparing the multi-fasta
- Steps to get all the IS110 and IS1111
    - one script frm the multifasta craetion, the rest is bash one liners: update_missing_sequences_fasta_headers.pl

### BLAST and filter    
blast_IS_complete.pbs 
filter_first_round_blast.pl

### Extract insertion site flanking sequence
extract_flanks_parallel.pl 
extract_flanks_make_chunks.pl  
extract_flanks_parallel_run.sh 
extract_flanks_concatenate.sh
extract_shorter_flanks.pl --> not needed, adjust RS script to take from centre of seq of any length and have params for size

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


