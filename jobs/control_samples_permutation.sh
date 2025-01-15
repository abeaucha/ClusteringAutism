#!/bin/bash
#SBATCH --job-name=control_samples_permutation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/control_samples_permutation_%j.out

# Activate virtual environment
source activate_venv_hpc.sh

control_samples_permutation.py