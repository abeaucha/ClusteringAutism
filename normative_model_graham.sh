#!/bin/bash
#SBATCH --job-name=normative_growth
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/normative_growth_%j.out

source activate_venv.sh

ti=$(date +"%T")
echo "Start time: $ti"
Rscript normative_growth_normalization.R \
	--imgdir /scratch/abeaucha/data/human/derivatives/POND_SickKids/jacobians/resolution_0.8/absolute/ \
	--demographics data/human/registration/DBM_input_demo_passedqc_wfile.csv \
	--outdir /scratch/abeaucha/data/human/derivatives/POND_SickKids/effect_sizes/698/resolution_0.8/absolute/ \
	--mask data/human/registration/reference_files/mask_0.8mm_0.1.mnc \
	--batch Site-Scanner \
	--nproc $SLURM_CPUS_PER_TASK
tf=$(date +"%T")
echo "End time: $tf"
