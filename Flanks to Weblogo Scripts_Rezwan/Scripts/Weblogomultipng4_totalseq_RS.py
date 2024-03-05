import os
from Bio import SeqIO
from weblogo import SeqList, LogoData, LogoOptions, LogoFormat, png_formatter, unambiguous_dna_alphabet

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
process_directory(
    r"C:\Users\Rezwan Siddiquee\Dropbox (Sydney Uni)\Transposases_RS\Target search Bioinformatics\Weblogo Generation\Weblogo Input_RS\20bp_flanks",
    r"C:\Users\Rezwan Siddiquee\Dropbox (Sydney Uni)\Transposases_RS\Target search Bioinformatics\Weblogo Generation\Weblogo Output_RS\20bp_Targets_RS"
)
