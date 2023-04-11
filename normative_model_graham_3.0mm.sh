#!/bin/bash
#SBATCH --job-name=normative_growth
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/normative_growth_3.0mm_%j.out

source activate_venv.sh

outdir=data/human/test/

ti=$(date +"%T")
echo "Start time: $ti"
#Rscript normative_growth_normalization.R \
#	--imgdir data/human/registration/v1/jacobians_resampled/resolution_3.0/absolute/ \
#	--demographics data/human/registration/v1/DBM_input_demo_passedqc_wfile.csv \
#	--outdir $SLURM_TMPDIR \
#	--mask data/human/registration/v1/reference_files/mask_3.0mm.mnc \
#	--batch Site-Scanner \
#	--nproc $SLURM_CPUS_PER_TASK
python3 normative_growth_test.py
tf=$(date +"%T")
echo "End time: $tf"

#cp -r $SLURM_TMPDIR $outdir
