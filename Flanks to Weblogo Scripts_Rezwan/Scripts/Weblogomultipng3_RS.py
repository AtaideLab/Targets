import os
from Bio import SeqIO
from weblogo import SeqList, LogoData, LogoOptions, LogoFormat, png_formatter, unambiguous_dna_alphabet

def create_weblogo(fasta_file, output_file, start=41, end=80):
    # Read sequences from fasta file
    seq_records = [record for record in SeqIO.parse(fasta_file, "fasta")]

    # Check if the sequence list is empty
    if not seq_records:
        print(f"Warning: No sequences found in {fasta_file}")
        return

    # Slice each sequence to include only positions 41 to 80 (Python uses 0-based indexing)
    sliced_sequences = [str(record.seq[start-1:end]) for record in seq_records]

    # Create SeqList from sliced_sequences
    sequences = SeqList(sliced_sequences, unambiguous_dna_alphabet)

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
process_directory(
    r"Place input folder directory here",
    r"Place output folder directory here"
)
