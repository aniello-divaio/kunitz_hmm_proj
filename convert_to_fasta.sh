#!/bin/bash

# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "ERROR:"
    echo -e "This program takes one argument from the command line: input file name.\nExiting..."
    exit 1
fi

input_file="$1"
output_file="kunitz_sequences.fasta"

# Clear or create output file
> "$output_file"

# Process the file
awk -F',' '{
    # Remove quotes from all fields
    for (i = 1; i <= NF; i++) {
        gsub(/"/, "", $i);
    }

    # Save ID, sequence, and chain if present
    if ($1 != "") last_id = $1;
    if ($4 != "") last_seq = $4;
    if ($5 != "") last_chain = $5;

    # If PF00014 is in the annotation field
    if ($6 ~ /PF00014/) {
        header = ">" last_id;
        if (last_chain != "") {
            header = header "_" last_chain;
        }
        print header >> "'$output_file'";
        print last_seq >> "'$output_file'";
    }
}' "$input_file"

# Notify the user about the output file
echo "Output file generated: $output_file"
