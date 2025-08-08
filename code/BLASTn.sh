#!/bin/bash
#PBS -l select=1:ncpus=16:mem=8gb
#PBS -l walltime=48:00:00

#(how many fasta files do I have?)
grep -c ">" "../data/BLAST/query_seqs/farmkit_ESV_sequences.fasta"

# Move to the directory where the job was submitted
cd $PBS_O_WORKDIR

# Create an output directory if it doesn't exist
mkdir -p output

# Load necessary modules (if any)
module load blast+/2.11.1

# List of input files in the input_files directory
input_files=("input_files/farmkit_ESV_sequences.part_001.fasta" \
"input_files/farmkit_ESV_sequences.part_002.fasta" \
"input_files/farmkit_ESV_sequences.part_003.fasta" \
"input_files/farmkit_ESV_sequences.part_004.fasta" \
"input_files/farmkit_ESV_sequences.part_005.fasta" \
"input_files/farmkit_ESV_sequences.part_006.fasta" \
"input_files/farmkit_ESV_sequences.part_007.fasta")

# Get the array task index (PBS_ARRAYID) to select the corresponding input file
input_file="${input_files[$PBS_ARRAYID-1]}"

# Run blastn for the corresponding file and save the result in the output directory
blastn -query "$input_file" -db /rds/general/user/hjt24/home/blast_db/silva_db/silva_db -out "output/$(basename "$input_file").results.txt" -outfmt 6 -num_threads 16
