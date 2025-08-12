#!/bin/bash

csv_file="orthogene_9species.csv"
output_dir="merge_7species"

# creat output folder
mkdir -p "$output_dir"

# read csv
tail -n +2 "$csv_file" | while IFS=',' read -r human_id cow_id pig_id mouse_id crabmacaque_id sheep_id horse_id rat_id marmoset_id gene_symbol; do
    # clean gene_symbol 
    gene_symbol=$(echo "$gene_symbol" | tr -d '\r' | tr -d '\n' | tr -d '[:space:]')

    output_file="${output_dir}/${gene_symbol}.fa"
    
    all_species_present=true

    # Extract gene sequence
    for species in human mouse rat cow pig marmoset crab_macaque; do
        case $species in
            human)
                gene_id=$human_id
                ;;
            mouse)
                gene_id=$mouse_id
                ;;
            rat)
                gene_id=$rat_id
                ;;
            cow)
                gene_id=$cow_id
                ;;
            pig)
                gene_id=$pig_id
                ;;
            marmoset)
                gene_id=$marmoset_id
                ;;
            crab_macaque)
                gene_id=$crabmacaque_id
                ;;
        esac
        
        file_path="/home/yf358/ref_backup/test_v2/${species}/genes_output/${gene_id}.fasta"
        if [ -f "$file_path" ]; then
            echo "Processing $file_path"
            cat "$file_path" >> "$output_file"
        else
            echo "Warning: File $file_path not found."
            all_species_present=false
        fi
    done

    if [ "$all_species_present" = false ]; then
        rm -f "$output_file"
    fi
done < <(cat "$csv_file"; echo)

echo "All files have been processed and merged into the 'merge' directory."