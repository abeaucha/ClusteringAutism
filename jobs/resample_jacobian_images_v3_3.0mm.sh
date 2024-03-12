#!/bin/bash
#SBATCH --job-name=resample_jacobian_images_v3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/resample_jacobian_images_v3_3.0mm_%j.out

# This pipeline ran in ___ minutes with --time=__:__:__
# Should be able to run it with a shorter walltime.

# Activate virtual environment
source activate_venv_hpc.sh

# Resample absolute images
resample_images.py \
--imgdir data/human/registration/v3/jacobians/absolute/smooth_minc/ \
--outdir data/human/registration/v3/jacobians_resampled/resolution_3.0/absolute/ \
--isostep 3.0 \
--nproc $SLURM_CPUS_PER_TASK

# Resample relative images
resample_images.py \
--imgdir data/human/registration/v3/jacobians/relative/smooth_minc/ \
--outdir data/human/registration/v3/jacobians_resampled/resolution_3.0/relative/ \
--isostep 3.0 \
--nproc $SLURM_CPUS_PER_TASK

deactivate
