#!/bin/bash
#SBATCH --job-name=metric_uncertainty
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=00:30:00
#SBATCH --chdir=/scratch/abeaucha/ClusteringAutism/main/figures/v3/

# Activate virtual environment
source activate_venv.sh

Rscript evaluate_cluster_metric_uncertainty.R

deactivate