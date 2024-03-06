#!/usr/bin/env python

import os
import sys
from Bio import SeqIO
from weblogo import SeqList, LogoData, LogoOptions, LogoFormat, png_formatter, unambiguous_dna_alphabet


# Input directory iwth fasta required argument:
if len(sys.argv) > 1:
    input_dir = sys.argv[1] 
else:
    print("No input fasta directory provided.\nPlease provide fasta directory as command-line argument to script.\nExiting.")
    sys.exit(1)



# Get intermediate and last parts of the input directory path
# Assumes the input directory follows the format: Output/Flanking_fastas_{filter_name}/{flank_size}
path_parts = os.path.normpath(input_dir).split(os.sep)
flank_size = path_parts[-1]
filter_name = path_parts[-2].split('_')[-2:]
output_dir_name = '_'.join(filter_name + [flank_size])

# Create the full path for the output directory
output_dir = os.path.join('Output', 'WebLogos', output_dir_name)
os.makedirs(output_dir, exist_ok=True)

print(f"Creating WebLogos on fastas in {input_dir}")
print(f"Writing WebLogos to {output_dir}")


def create_weblogo(fasta_file, output_file):
    # Read sequences from fasta file
    seq_records = [record for record in SeqIO.parse(fasta_file, "fasta")]

    # Check if the sequence list is empty
    if not seq_records:
        print(f"Warning: No sequences found in {fasta_file}")
        return

    # Create SeqList from seq_records
    sequences = SeqList([str(record.seq) for record in seq_records], unambiguous_dna_alphabet)

    # Create a multiple sequence alignment
    data = LogoData.from_seqs(sequences)

    # Generate options
    options = LogoOptions()
    options.title = "WebLogo Example"
    
    # Generate the logo
    format = LogoFormat(data, options)
    with open(output_file, "wb") as f:
        f.write(png_formatter(data, format))

def process_directory(input_dir, output_dir):
    # Loop over all files in the input directory
    for filename in os.listdir(input_dir):
        # Check if the file is a FASTA file
        if filename.endswith(".fasta"):
            # Construct the full paths for the input and output files
            fasta_file = os.path.join(input_dir, filename)
            png_file = os.path.join(output_dir, os.path.splitext(filename)[0] + ".png")
            
            print(f"Processing: {filename}")  # Print the file being processed
            create_weblogo(fasta_file, png_file)

# Usage
process_directory(input_dir, output_dir)
###
