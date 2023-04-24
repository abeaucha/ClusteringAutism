#!/bin/bash
#SBATCH --job-name=permute_cluster_similarity_v2_rel_1
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/permute_cluster_similarity_v2_0.8mm_rel_1_%j.out

source activate_venv.sh

#Can only do 4 permutations in 48 hours
ti=$(date)
echo "Start time: $ti"
python3 permute_cluster_similarity.py \
	--pipeline-dir data/cross_species/v2/ \
	--params-id 405 \
	--human-pipeline-dir data/human/derivatives/v2/ \
	--mouse-pipeline-dir data/mouse/derivatives/v2/ \
	--mouse-expr-dir data/mouse/expression/ \
	--human-expr-dir data/human/expression/ \
	--mouse-mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	--human-mask data/human/registration/v2/reference_files/mask_0.8mm.mnc \
	--human-microarray-coords data/human/expression/AHBA_microarray_coordinates_study_v2.csv \
	--npermutations 4 \
	--permutations-start 5 \
	--jacobians relative \
	--nproc $SLURM_CPUS_PER_TASK
tf=$(date)
echo "End time: $tf"