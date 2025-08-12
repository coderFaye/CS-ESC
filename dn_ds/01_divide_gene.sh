#!/bin/bash

# Input FASTA file
input_fasta="Homo_sapiens.GRCh38.cds.test.fa"

# Create a directory to store output files
mkdir -p genes_output

# Variables to keep track of the longest transcript
current_gene=""
current_gene_file=""
longest_seq=""
longest_length=0

# Function to write the sequence to the file
write_sequence() {
    if [[ -n "$current_gene" && -n "$longest_seq" ]]; then
        # Remove the decimal point part of the Ensembl gene ID before writing
        gene_id=$(echo "$current_gene" | cut -d'.' -f1)
        # Output the longest sequence to the gene's file
        echo ">human" > "genes_output/${gene_id}.fasta"
        echo "$longest_seq" >> "genes_output/${gene_id}.fasta"
    fi
}

# Initialize variables
sequence=""
gene_id=""
gene_name=""

# Read the FASTA file
while read -r line; do
    if [[ $line == ">"* ]]; then
        # We reached a new transcript line, so process the previous one

        # If this is not the first gene, check if we should write the longest transcript to the file
        if [[ -n "$current_gene" ]]; then
            if [[ ${#sequence} -gt $longest_length ]]; then
                longest_seq="$sequence"
                longest_length=${#sequence}
            fi
        fi

        # Extract the Ensembl gene ID from the current line
        gene_id=$(echo "$line" | grep -oP "gene:\K[^ ]+")
        gene_name=$(echo "$gene_id" | cut -d'.' -f1)

        # If this is a new gene, write the previous longest transcript to the file
        if [[ "$gene_name" != "$current_gene" ]]; then
            write_sequence
            # Reset for the new gene
            current_gene="$gene_name"
            longest_seq=""
            longest_length=0
        fi

        # Reset the sequence for the new gene
        sequence=""
    else
        # Collect the sequence lines
        sequence+="$line"
    fi
done < "$input_fasta"

# Write the final gene's longest sequence
write_sequence

echo "Separation of genes into separate FASTA files is complete."