#!/bin/bash
#SBATCH --job-name=permute_cluster_similarity
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/permute_cluster_similarity_v1_0.8mm_%j.out

source activate_venv.sh

ti=$(date +"%T")
echo "Start time: $ti"
python3 permute_cluster_similarity.py \
	--pipeline-dir data/cross_species/v1/ \
	--params-id 405 \
	--human-pipeline-dir data/human/derivatives/v1/ \
	--mouse-pipeline-dir data/mouse/derivatives/v2/ \
	--mouse-expr-dir data/mouse/expression/ \
	--human-expr-dir data/human/expression/ \
	--mouse-mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	--human-mask data/human/registration/v1/reference_files/mask_0.8mm.mnc \
	--human-microarray-coords data/human/expression/AHBA_microarray_coordinates_study_v1.csv \
	--npermutations 10 \
	--nproc $SLURM_CPUS_PER_TASK
tf=$(date +"%T")
echo "End time: $tf"