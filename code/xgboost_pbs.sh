#!/bin/bash
#PBS -N xgboost_316k
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=16:mem=128gb

module load R/4.4.2-gfbf-2024a

cd $PBS_O_WORKDIR

# Force export cores
export NCPUS=16
export OMP_NUM_THREADS=16

echo "Starting XGBoost analysis with $NCPUS cores"
Rscript xgboost_hpc.R