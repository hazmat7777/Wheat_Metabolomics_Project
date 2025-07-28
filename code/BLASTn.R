#!/bin/bash
#PBS -N BLAST_search              # Job name
#PBS -o blast_output.txt          # Standard output file
#PBS -e blast_error.txt           # Standard error file
#PBS -l walltime=72:00:00         # Max walltime (72 hours)
#PBS -l nodes=1:ppn=16            # Request 1 node with 16 cores
#PBS -l mem=64gb                  # Memory request per node
#PBS -q medium72                  # Queue partition for medium jobs
#PBS -t 1-7                       # Job array for 7 input files

# Move to the directory where the job was submitted
cd $PBS_O_WORKDIR

# Create an output directory if it doesn't exist
mkdir -p output

# Load necessary modules (if any)
module load ncbi-blast/2.17.0

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
blastn -query "$input_file" -db /path/to/silva_db -out "output/$(basename "$input_file").results.txt" -outfmt 6 -num_threads 16
