#!/bin/bash
#SBATCH --job-name=normative_growth
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/normative_growth_0.8mm_%j.out

source activate_venv.sh

ti=$(date +"%T")
echo "Start time: $ti"
python3 normative_growth_test.py
tf=$(date +"%T")
echo "End time: $tf"
