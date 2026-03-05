#!/bin/bash
#SBATCH --job-name=metric_uncertainty
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=12:00:00
#SBATCH --chdir=/scratch/abeaucha/ClusteringAutism/main/

# Activate virtual environment
source activate_venv.sh

Rscript figures/v3/evaluate_cluster_metric_uncertainty.R

deactivate
