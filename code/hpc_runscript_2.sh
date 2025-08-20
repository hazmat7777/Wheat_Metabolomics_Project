#!/bin/bash
#PBS -N rf2_316k_taxa
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=16:mem=128gb

# Load required modules
module load R/4.4.2-gfbf-2024a

# Change to working directory
cd $PBS_O_WORKDIR

# Print job info
echo "Job started at $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $PBS_JOBID"

# CRITICAL: Force export the number of CPUs
# PBS might use different variable names on different systems
export NCPUS=16
export OMP_NUM_THREADS=16
export PBS_NP=16

echo "Forcing NCPUS to: $NCPUS"

# Run the R script
echo "Starting Random Forest analysis on 316k taxa..."
Rscript microb_rf_hpc_2.R

echo "Job completed at $(date)"