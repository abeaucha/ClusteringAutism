#!/usr/bin/env bash
#SBATCH --job-name=slurm_test
#SBATCH -N 1 
#SBATCH --cpus-per-task 8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/%j.out

export QBATCH_SYSTEM=slurm

source activate_hpfenv.sh

echo $(which R)
echo $(which python3)
echo $MODULEPATH
echo $PYTHONPATH

# Rscript slurm_test.R
python3 slurm_dev.py
