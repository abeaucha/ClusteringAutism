#!/bin/bash
#SBATCH --job-name=compute_cluster_similarity
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/permute_cluster_similarity_%j.out

source activate_venv.sh

python3 permute_cluster_similarity.py \
  --human-pipeline-dir data/human/derivatives/v1/ \
  --mouse-pipeline-dir data/mouse/derivatives/ \
  --mouse-expr-dir data/mouse/expression/ \
	--human-expr-dir data/human/expression/ \
	--mouse-mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	--human-mask data/human/registration/v1/reference_files/mask_3.0mm.mnc \
	--human-microarray-coords data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
	--output-dir data/cross_species/v1/ \
	--human-cluster-params-id 488-834 \
	--human-dataset POND_SickKids \
	--mouse-dataset Models_135 \
	--npermutations 50 \
	--cluster-map-method mean \
	--sim-gene-space average-latent-space \
	--sim-n-latent-spacs 50 \
	--nproc $SLURM_CPUS_PER_TASK