#!/bin/bash
#SBATCH --job-name=normative_growth
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --tmp=30G
#SBATCH --time=6:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/normative_growth_%j.out

source activate_venv.sh

ti=$(date +"%T")
echo "Start time: $ti"
Rscript normative_growth_normalization_test.R \
	--imgdir data/human/registration/v1/jacobians_resampled/resolution_0.8/absolute/ \
	--demographics data/human/registration/v1/DBM_input_demo_passedqc_wfile.csv \
	--outdir $SLURM_TMPDIR \
	--mask data/human/registration/v1/reference_files/mask_0.8mm.mnc \
	--batch Site-Scanner \
	--nproc $SLURM_CPUS_PER_TASK
tf=$(date +"%T")
echo "End time: $tf"
