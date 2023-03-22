#!/usr/bin/env bash
#SBATCH -N 1 -c 8
#SBATCH --mem=32G
#SBATCH --time=6:00:00

export QBATCH_SYSTEM=slurm

cd /hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main

source activate_hpfenv.sh

echo $MODULEPATH

python3 slurm_dev.py