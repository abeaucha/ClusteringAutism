#!/bin/bash
#SBATCH --job-name=compute_cluster_similarity
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/compute_cluster_similarity_%j.out

source activate_venv.sh

outdir=data/cross_species/v1/Models_135-POND_SickKids/similarity/

python3 compute_cluster_similarity.py \
	--pipeline-dir $outdir \
	--mouse-cluster-dir data/mouse/derivatives/Models_135/cluster_maps/956/ \
	--human-cluster-dir data/human/derivatives/v1/POND_SickKids/cluster_maps/488-834-817/ \
	--mouse-expr-dir data/mouse/expression/ \
	--human-expr-dir data/human/expression/ \
	--mouse-mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	--human-mask data/human/registration/v1/reference_files/mask_3.0mm.mnc \
	--mouse-resolution 0.2 \
	--human-resolution 3.0 \
	--human-microarray-coords data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
	--gene-space average-latent-space \
	--n-latent-spaces 5 \
	--nproc $SLURM_CPUS_PER_TASK
