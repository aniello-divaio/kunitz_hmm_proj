"""
get_seq.py

This script extracts specific protein sequences from a UniProt-formatted FASTA file, 
based on a list of UniProt accession IDs provided in a separate text file.

The output is a FASTA file named like the ID file, but with '.fasta' extension.

Usage:
    python get_seq.py <id_list.txt> <input_fasta.fasta>
"""

from sys import argv
from Bio import SeqIO
import os

def get_ids(file):
    with open(file, "r") as f:
        return {line.strip() for line in f if line.strip()}

def get_seq(ids_to_get, input_fasta, output_fasta):
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            try:
                uniprot_id = record.id.split('|')[1] if '|' in record.id else record.id
            except IndexError:
                continue
            if uniprot_id in ids_to_get:
                SeqIO.write(record, outfile, "fasta")

if __name__ == "__main__":
    if len(argv) != 3:
        print("ERROR: This program takes an ID file and a FASTA file as arguments.")
        print("Usage: python get_seq.py <id_list.txt> <input_fasta.fasta>")
        exit(1)
    else:
        id_file = argv[1]
        input_fasta = argv[2]

        # Output will be the same as the ID file but with '.fasta' extension
        base_name = os.path.splitext(id_file)[0]
        output_fasta = f"{base_name}.fasta"

        ids = get_ids(id_file)
        get_seq(ids, input_fasta, output_fasta)

        print(f"Written output to: {output_fasta}")
