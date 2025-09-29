#!/bin/bash
#SBATCH --job-name=compute_cluster_similarity
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main/
#SBATCH --output=logs/compute_cluster_similarity_PONDSK_%j.out
##SBATCH --dependency=afterok:

# Activate virtual environment
source activate_venv.sh

# Pipeline registry directory
REGISTRY="compute_cluster_similarity_registry_${SLURM_JOB_ID}"

# Execute pipeline
# Maximum of 300 jobs allowed on HPF (??)
compute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--species human mouse \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--input-params-ids 700 107 \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--gene-space vae-latent-space \
--jacobians absolute relative \
--execution slurm \
--nproc 600 \
--slurm-mem 16G \
--slurm-time 8:00:00

deactivate