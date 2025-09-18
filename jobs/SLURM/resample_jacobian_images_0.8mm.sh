#!/bin/bash
#SBATCH --job-name=resample_jacobian_images_v3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/resample_jacobian_images_v3_0.8mm_%j.out

# This pipeline ran in 40 minutes with --time=12:00:00
# Should be able to run it with a shorter walltime.

# Activate virtual environment
source activate_venv_hpc.sh

# Resample absolute images
resample_images.py \
--imgdir data/human/registration/v3/jacobians/absolute/smooth_minc/ \
--outdir data/human/registration/v3/jacobians_resampled/resolution_0.8/absolute/ \
--isostep 0.8 \
--nproc $SLURM_CPUS_PER_TASK

# Resample relative images
resample_images.py \
--imgdir data/human/registration/v3/jacobians/relative/smooth_minc/ \
--outdir data/human/registration/v3/jacobians_resampled/resolution_0.8/relative/ \
--isostep 0.8 \
--nproc $SLURM_CPUS_PER_TASK

deactivate
