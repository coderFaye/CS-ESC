########translatorx - pal2nal - oneline_formating - phylip_formating - raxml - codeml########

### translatorx (alignment)
chmod +x /home/yf358/ref_backup/test/perl/1_translator_x/translatorX.pl
perl /home/yf358/ref_backup/test/perl/1_translator_x/translatorX.pl \
-i /home/yf358/ref_backup/test_v2/orthology/merge_7species/ZZZ3.fa \
-o /home/yf358/ref_backup/test_v2/translatorx/ZZZ3.out \
-p F

# Loop over all .fa files in the folder (using packaged script: translatorx_parelle.sh)
mkdir /home/yf358/ref_backup/test_v2/translatorx_logs
mkdir /home/yf358/ref_backup/test_v2/translatorx_p
chmod +x /home/yf358/ref_backup/test_v2/translatorx_parelle.sh
nohup sh /home/yf358/ref_backup/test_v2/translatorx_parelle.sh &

# Number of parallel groups
num_groups=20
# Create an array of all gene files
gene_files=(/home/yf358/ref_backup/test_v2/orthology/merge_7species/*.fa)
# Calculate the number of genes per group
num_genes_per_group=$(( (${#gene_files[@]} + num_groups - 1) / num_groups ))
# Loop over the number of groups
for (( group=0; group<num_groups; group++ )); do
    # Calculate the start and end index for this group
    start_index=$(( group * num_genes_per_group ))
    end_index=$(( (group + 1) * num_genes_per_group - 1 ))

    # Make sure the end index doesn't exceed the total number of genes
    if [[ $end_index -ge ${#gene_files[@]} ]]; then
        end_index=$((${#gene_files[@]} - 1))
    fi

    # Create a subset of gene files for this group
    group_files=("${gene_files[@]:$start_index:$((end_index - start_index + 1))}")

    # Create a temporary script to run in the background
    temp_script="/tmp/run_group_$group.sh"
    echo "#!/bin/bash" > "$temp_script"
    
    for gene_file in "${group_files[@]}"; do
        # Extract the gene name from the file path
        gene_name=$(basename "$gene_file" .fa)

        # Append the command to the temp script
        echo "perl /home/yf358/ref_backup/test/perl/1_translator_x/translatorX.pl \
        -i '$gene_file' \
        -o '/home/yf358/ref_backup/test_v2/translatorx_p/${gene_name}.out' \
        -p F" >> "$temp_script"

        echo "echo 'Processed: $gene_name'" >> "$temp_script"
    done

    # Run the temp script in the background with nohup
    nohup bash "$temp_script" > "/home/yf358/ref_backup/test_v2/translatorx_logs/group_$group.log" 2>&1 &

    echo "Started group $group, processing files from $start_index to $end_index"
done

@parameter:
# -p F \ use MAFFT
@output_file:
# ZZZ3.out.nt_ali.fasta #nuecleotide alignment
# ZZZ3.out.aa_ali.fasta #amino acid alignment



### pal2nal (convert the alignments in codon data to double check the nucleotide alignment is the same)
chmod +x /home/yf358/ref_backup/test/perl/2_pal2nal/pal2nal.pl
/home/yf358/ref_backup/test/perl/2_pal2nal/pal2nal.pl \
/home/yf358/ref_backup/test_v2/translatorx/ZZZ3.out.aa_ali.fasta /home/yf358/ref_backup/test_v2/translatorx/ZZZ3.out.nt_ali.fasta \
-output fasta > /home/yf358/ref_backup/test_v2/pal2nal/ZZZ3_pal2nal_mafft_out.fasta

# Loop over all .out.aa_ali.fasta files in the translatorx output folder (using packaged script pal2nal_paralle.sh)
mkdir /home/yf358/ref_backup/test_v2/pal2nal_logs
chmod +x /home/yf358/ref_backup/test_v2/pal2nal_paralle.sh
nohup sh /home/yf358/ref_backup/test_v2/pal2nal_paralle.sh &


# Number of parallel groups
num_groups=20
# Create an array of all .out.aa_ali.fasta files in the translatorx output folder
aa_files=(/home/yf358/ref_backup/test_v2/translatorx_p/*.out.aa_ali.fasta)
# Calculate the number of files per group
num_files_per_group=$(( (${#aa_files[@]} + num_groups - 1) / num_groups ))
# Loop over the number of groups
for (( group=0; group<num_groups; group++ )); do
    # Calculate the start and end index for this group
    start_index=$(( group * num_files_per_group ))
    end_index=$(( (group + 1) * num_files_per_group - 1 ))

    # Make sure the end index doesn't exceed the total number of files
    if [[ $end_index -ge ${#aa_files[@]} ]]; then
        end_index=$((${#aa_files[@]} - 1))
    fi

    # Create a subset of aa_files for this group
    group_files=("${aa_files[@]:$start_index:$((end_index - start_index + 1))}")

    # Create a temporary script to run in the background
    temp_script="/tmp/run_pal2nal_group_$group.sh"
    echo "#!/bin/bash" > "$temp_script"
    
    for aa_file in "${group_files[@]}"; do
        # Extract the gene name from the file path
        gene_name=$(basename "$aa_file" .out.aa_ali.fasta)

        # Define the corresponding nucleotide alignment file
        nt_file="/home/yf358/ref_backup/test_v2/translatorx_p/${gene_name}.out.nt_ali.fasta"

        # Define the output file for pal2nal
        output_file="/home/yf358/ref_backup/test_v2/pal2nal/${gene_name}_pal2nal_mafft_out.fasta"

        # Append the pal2nal command to the temp script
        echo "/home/yf358/ref_backup/test/perl/2_pal2nal/pal2nal.pl \
        '$aa_file' \
        '$nt_file' \
        -output fasta > '$output_file'" >> "$temp_script"

        echo "echo 'Processed PAL2NAL for: $gene_name'" >> "$temp_script"
    done

    # Run the temp script in the background with nohup
    nohup bash "$temp_script" > "/home/yf358/ref_backup/test_v2/pal2nal_logs/group_$group.log" 2>&1 &

    echo "Started group $group, processing files from $start_index to $end_index"
done



### Convert output file to one line fasta file
chmod +x /home/yf358/ref_backup/test/perl/2_pal2nal/one_line_fasta.pl
/home/yf358/ref_backup/test/perl/one_line_fasta.pl ZZZ3_pal2nal_mafft_out.fasta

# loop for all the files (using change_to_one_line_parallel.sh)
mkdir /home/yf358/ref_backup/test_v2/one_line_fasta_logs
chmod +x /home/yf358/ref_backup/test_v2/change_to_one_line_parallel.sh
nohup sh /home/yf358/ref_backup/test_v2/change_to_one_line_parallel.sh &

# Number of parallel groups
num_groups=20
# Create an array of all _pal2nal_mafft_out.fasta files in the pal2nal output folder
pal2nal_files=(/home/yf358/ref_backup/test_v2/pal2nal/*_pal2nal_mafft_out.fasta)
# Calculate the number of files per group
num_files_per_group=$(( (${#pal2nal_files[@]} + num_groups - 1) / num_groups ))
# Loop over the number of groups
for (( group=0; group<num_groups; group++ )); do
    # Calculate the start and end index for this group
    start_index=$(( group * num_files_per_group ))
    end_index=$(( (group + 1) * num_files_per_group - 1 ))

    # Make sure the end index doesn't exceed the total number of files
    if [[ $end_index -ge ${#pal2nal_files[@]} ]]; then
        end_index=$((${#pal2nal_files[@]} - 1))
    fi

    # Create a subset of pal2nal_files for this group
    group_files=("${pal2nal_files[@]:$start_index:$((end_index - start_index + 1))}")

    # Create a temporary script to run in the background
    temp_script="/tmp/run_one_line_fasta_group_$group.sh"
    echo "#!/bin/bash" > "$temp_script"
    
    for pal2nal_file in "${group_files[@]}"; do
        # Extract the gene name from the file path
        gene_name=$(basename "$pal2nal_file" _pal2nal_mafft_out.fasta)

        # Append the one_line_fasta.pl command to the temp script
        echo "/home/yf358/ref_backup/test/perl/2_pal2nal/one_line_fasta.pl '$pal2nal_file'" >> "$temp_script"

        echo "echo 'Processed one_line_fasta for: $gene_name'" >> "$temp_script"
    done

    # Run the temp script in the background with nohup
    nohup bash "$temp_script" > "/home/yf358/ref_backup/test_v2/one_line_fasta_logs/group_$group.log" 2>&1 &

    echo "Started group $group, processing files from $start_index to $end_index"
done



### convert fasta in phylip
chmod +x /home/yf358/ref_backup/test/perl/2_pal2nal/FASTAtoPHYL.pl
num=$( grep '>' ZZZ3_pal2nal_mafft_out_one_line.fa | wc -l )
len=$( sed -n '2,2p' ZZZ3_pal2nal_mafft_out_one_line.fa | sed 's/\r//' | sed 's/\n//' | wc -L )
perl /home/yf358/ref_backup/test/perl/2_pal2nal/FASTAtoPHYL.pl ZZZ3_pal2nal_mafft_out_one_line.fa $num $len 
mv ZZZ3_pal2nal_mafft_out_one_line.fa.phy /home/yf358/ref_backup/test_v2/phylip_output/ZZZ3_mafft.phy

# loop for all the genes (using to_phy_parallel.sh)
mkdir /home/yf358/ref_backup/test_v2/phylip_logs
mkdir /home/yf358/ref_backup/test_v2/phylip_output
chmod +x /home/yf358/ref_backup/test_v2/to_phy_parallel.sh
cd /home/yf358/ref_backup/test_v2/phylip_output
nohup sh /home/yf358/ref_backup/test_v2/to_phy_parallel.sh &

# Number of parallel groups
num_groups=20
# Define the output directory for .phy files
output_dir="/home/yf358/ref_backup/test_v2/phylip_output"
# Create an array of all _pal2nal_mafft_out_one_line.fa files in the pal2nal folder
fasta_files=(/home/yf358/ref_backup/test_v2/pal2nal/*_pal2nal_mafft_out_one_line.fa)
# Calculate the number of files per group
num_files_per_group=$(( (${#fasta_files[@]} + num_groups - 1) / num_groups ))
# Loop over the number of groups
for (( group=0; group<num_groups; group++ )); do
    # Calculate the start and end index for this group
    start_index=$(( group * num_files_per_group ))
    end_index=$(( (group + 1) * num_files_per_group - 1 ))

    # Make sure the end index doesn't exceed the total number of files
    if [[ $end_index -ge ${#fasta_files[@]} ]]; then
        end_index=$((${#fasta_files[@]} - 1))
    fi

    # Create a subset of fasta_files for this group
    group_files=("${fasta_files[@]:$start_index:$((end_index - start_index + 1))}")

    # Create a temporary script to run in the background
    temp_script="/tmp/run_fasta_to_phylip_group_$group.sh"
    echo "#!/bin/bash" > "$temp_script"
    
    for fasta_file in "${group_files[@]}"; do
        # Extract the gene name by removing the _pal2nal_mafft_out_one_line.fa suffix
        gene_name=$(basename "$fasta_file" _pal2nal_mafft_out_one_line.fa)

        # Count the number of sequences (num) in the FASTA file
        num=$(grep '>' "$fasta_file" | wc -l)

        # Get the length of the sequences (len) by looking at the second line of the file
        len=$(sed -n '2,2p' "$fasta_file" | sed 's/\r//' | sed 's/\n//' | wc -L)

        # Append the FASTAtoPHYL.pl command to the temp script
        echo "perl /home/yf358/ref_backup/test/perl/2_pal2nal/FASTAtoPHYL.pl '$fasta_file' $num $len" >> "$temp_script"

        # Append the mv command to move the .phy file to the output directory
        echo "mv '${fasta_file}.phy' '$output_dir/${gene_name}_mafft.phy'" >> "$temp_script"

        echo "echo 'Processed PHYLIP conversion for: $gene_name'" >> "$temp_script"
    done

    # Run the temp script in the background with nohup
    nohup bash "$temp_script" > "/home/yf358/ref_backup/test_v2/phylip_logs/group_$group.log" 2>&1 &

    echo "Started group $group, processing files from $start_index to $end_index"
done

# loop to change file name
# Loop over all the .fa.phy files in the directory
for file in /home/yf358/ref_backup/test_v2/phylip_output/*_pal2nal_mafft_out_one_line.fa.phy; do
    # Extract the base gene name by removing the _pal2nal_mafft_out_one_line.fa.phy suffix
    gene_name=$(basename "$file" _pal2nal_mafft_out_one_line.fa.phy)

    # Define the new filename
    new_filename="/home/yf358/ref_backup/test_v2/phylip_output/${gene_name}.phy"

    # Rename the file
    mv "$file" "$new_filename"

    echo "Renamed: $file -> $new_filename"
done



### raxml (gene tree)
/programs/RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -p 12345 -# 100 -x 12345 -# 500 -s ZZZ3.phy -n ZZZ3
printf "7  1\n" > ./ZZZ3.tree
sed 's/\:[0-9]*\.[0-9]*//g' ./RAxML_bestTree.ZZZ3 >> ./ZZZ3.tree

# Loop through each .phy file in the folder (using)
cd /home/yf358/ref_backup/test_v2/phylip_output
mkdir /home/yf358/ref_backup/test_v2/raxml_logs
chmod +x /home/yf358/ref_backup/test_v2/raxml_paralle.sh
nohup sh /home/yf358/ref_backup/test_v2/raxml_paralle.sh &

# Number of parallel groups
num_groups=100
# Define the folder where the .phy files are located
phy_dir="/home/yf358/ref_backup/test_v2/phylip_output"
# Create an array of all .phy files in the folder
phy_files=("$phy_dir"/*.phy)
# Calculate the number of files per group
num_files_per_group=$(( (${#phy_files[@]} + num_groups - 1) / num_groups ))
# Loop over the number of groups
for (( group=0; group<num_groups; group++ )); do
    # Calculate the start and end index for this group
    start_index=$(( group * num_files_per_group ))
    end_index=$(( (group + 1) * num_files_per_group - 1 ))

    # Make sure the end index doesn't exceed the total number of files
    if [[ $end_index -ge ${#phy_files[@]} ]]; then
        end_index=$((${#phy_files[@]} - 1))
    fi

    # Create a subset of phy_files for this group
    group_files=("${phy_files[@]:$start_index:$((end_index - start_index + 1))}")

    # Create a temporary script to run in the background
    temp_script="/tmp/run_raxml_group_$group.sh"
    echo "#!/bin/bash" > "$temp_script"
    
    for file in "${group_files[@]}"; do
        # Extract the gene name by removing the .phy extension
        GENE_NAME=$(basename "$file" .phy)

        # Append the RAxML command to the temp script
        echo "/programs/RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 \
          -f a \
          -m GTRGAMMA \
          -p 12345 \
          -# 100 \
          -x 12345 \
          -# 500 \
          -s '$file' \
          -n '$GENE_NAME'" >> "$temp_script"

        # Append the commands to create the tree file and remove branch lengths
        echo "printf '7  1\n' > './${GENE_NAME}.tree'" >> "$temp_script"
        echo "sed 's/:[0-9]*\.[0-9]*//g' './RAxML_bestTree.${GENE_NAME}' >> './${GENE_NAME}.tree'" >> "$temp_script"
        echo "echo 'Processed RAxML for: $GENE_NAME'" >> "$temp_script"
    done

    # Run the temp script in the background with nohup
    nohup bash "$temp_script" > "/home/yf358/ref_backup/test_v2/raxml_logs/group_$group.log" 2>&1 &

    echo "Started group $group, processing files from $start_index to $end_index"
done
# kill the loop
ps aux | grep run_raxml_group
pkill -f run_raxml_group

@parameter:
# -f a  select algorithm:  "-f a": rapid Bootstrap analysis and search for best-scoring ML tree in one program run
# -m GTRGAMMA Model of Binary (Morphological), Nucleotide, Multi-State, or Amino Acid Substitution: "-m GTRGAMMA[X]"     : GTR + Optimization of substitution rates + GAMMA model of rate
#heterogeneity (alpha parameter will be estimated).
# -p Specify a random number seed for the parsimony inferences. This allows you to reproduce your results
#and will help me debug the program.
# -x Specify an integer number (random seed) and turn on rapid bootstrapping
# -s Specify the name of the alignment data file in PHYLIP format
# -o Specify the name of a single outgroup or a comma-separated list of outgroups, eg "-o Rat"
#or "-o Rat,Mouse", in case that multiple outgroups are not monophyletic the first name
#in the list will be selected as outgroup, don't leave spaces between taxon names!



### codeml
# prepare the control file for each gene
cd /home/yf358/ref_backup/test_v2/Model_M0
phy_dir="/home/yf358/ref_backup/test_v2/phylip_output"
# Loop over all .phy files in the directory
for phy_file in "$phy_dir"/*.phy; do
    # Extract the gene name (remove the directory path and .phy extension)
    gene_name=$(basename "$phy_file" .phy)
    
    # Copy the template control file to a new file named after the gene
    cp codeml-M0.ctl "${gene_name}.ctl"

    # Substitute the seqfile, treefile, and outfile in the new control file
    sed -i "s|../../Mx_aln.phy|$phy_dir/${gene_name}.phy|" "${gene_name}.ctl"
    sed -i "s|../../Mx_unroot.tree|$phy_dir/${gene_name}.tree|" "${gene_name}.ctl"
    sed -i "s|out_M0.txt|out_${gene_name}.txt|" "${gene_name}.ctl"

    echo "Prepared control file for gene: ${gene_name}"
done

# run the codeml for one gene
export PATH=/programs/paml-4.10.6/bin_force:$PATH
codeml ZZZ3.ctl | tee logfile_ZZZ3.txt
# Remove unnecessary files 
rm 2N*
# extract the estimated w in an output file
grep 'omega ' out_ZZZ3.txt | sed 's/..*= *//' >> omega_ZZZ3.txt 

# run the codeml for all the genes
export PATH=/programs/paml-4.10.6/bin_force:$PATH
mkdir /home/yf358/ref_backup/test_v2/codeml_logs
chmod +x /home/yf358/ref_backup/test_v2/Model_M0/codeml_paralle.sh
nohup sh /home/yf358/ref_backup/test_v2/Model_M0/codeml_paralle.sh &
chmod +x all_omega.sh
sh all_omega.sh

# Number of parallel groups
num_groups=100
# Define the folder where the .ctl files are located
ctl_dir="/home/yf358/ref_backup/test_v2/Model_M0"
# Create an array of all .ctl files in the folder
ctl_files=("$ctl_dir"/*.ctl)
# Calculate the number of files per group
num_files_per_group=$(( (${#ctl_files[@]} + num_groups - 1) / num_groups ))
# Loop over the number of groups
for (( group=0; group<num_groups; group++ )); do
    # Calculate the start and end index for this group
    start_index=$(( group * num_files_per_group ))
    end_index=$(( (group + 1) * num_files_per_group - 1 ))

    # Make sure the end index doesn't exceed the total number of files
    if [[ $end_index -ge ${#ctl_files[@]} ]]; then
        end_index=$((${#ctl_files[@]} - 1))
    fi

    # Create a subset of ctl_files for this group
    group_files=("${ctl_files[@]:$start_index:$((end_index - start_index + 1))}")

    # Create a temporary script to run in the background
    temp_script="/tmp/run_codeml_group_$group.sh"
    echo "#!/bin/bash" > "$temp_script"
    
    for ctl_file in "${group_files[@]}"; do
        # Extract the gene name by removing the .ctl extension
        gene_name=$(basename "$ctl_file" .ctl)

        # Append the codeml command to the temp script
        echo "codeml '$ctl_file' | tee logfile_${gene_name}.txt" >> "$temp_script"

        # Append the commands to remove unnecessary files and extract omega
        echo "rm 2N*" >> "$temp_script"
        echo "grep 'omega ' 'out_${gene_name}.txt' | sed 's/..*= *//' >> omega_${gene_name}.txt" >> "$temp_script"

        echo "echo 'Processed CODEML for: $gene_name'" >> "$temp_script"
    done

    # Run the temp script in the background with nohup
    nohup bash "$temp_script" > "/home/yf358/ref_backup/test_v2/codeml_logs/group_$group.log" 2>&1 &

    echo "Started group $group, processing files from $start_index to $end_index"
done