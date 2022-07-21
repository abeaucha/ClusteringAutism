#!/bin/bash

# ----------------------------------------------------------------------------
# pipeline_human.sh
# Author: Antoine Beauchamp
# 
# Pipeline to process and cluster human imaging data.
# Pipeline steps:
# 1. Extract jacobian images from compressed files
# 2. Calculate participant voxel-wise effect sizes
# 3. Downsample template, mask, and effect size images
# 4. Build effect sizes matrices from images
# 5. Cluster participants based on effect sizes
# 6. Generate representative cluster effect size maps


#Activate virtual environment 
source activate_venv.sh

#Option to use POND subset
dataset='ALL'
# dataset=POND

echo "Dataset is ${dataset}"

if [ -z $dataset ];
then
	human_dir=data/human/
else
	if [ $dataset = 'POND' ];
	then
		human_dir=data/human/POND/
	else
		human_dir=data/human/
	fi
fi


# Extract compressed jacobian images -----------------------------------------

#Extract compressed absolute jacobians and convert to minc
#python3 extract_zipped_images.py \
#	--imgdir data/human/registration/jacobians/absolute/smooth/ \
#	--parallel true \
#	--nproc 8

#Extract compressed relative jacobians and convert to minc
#python3 extract_zipped_images.py \
#	--imgdir data/human/registration/jacobians/relative/smooth/ \
#	--parallel true \
#	--nproc 8


# Calculate effect sizes -----------------------------------------------------

#Calculate human effect sizes based using absolute volumes
echo "Calculating human effect sizes using absolute jacobians..."
Rscript calculate_human_effect_sizes.R \
	--demographics data/human/registration/DBM_input_demo_passedqc.csv \
	--imgdir data/human/registration/jacobians/absolute/smooth/ \
	--outdir ${human_dir}/effect_sizes/absolute/resolution_0.5/ \
	--maskfile data/human/registration/reference_files/mask.mnc \
	--ncontrols 10 \
	--threshold 5 \
	--dataset $dataset \
	--parallel true \
	--nproc 4

#Calculate human effect sizes based using relative volumes
echo "Calculating human effect sizes using relative jacobians..."
Rscript calculate_human_effect_sizes.R \
	--demographics data/human/registration/DBM_input_demo_passedqc.csv \
	--imgdir data/human/registration/jacobians/relative/smooth/ \
	--outdir ${human_dir}/effect_sizes/relative/resolution_0.5/ \
	--maskfile data/human/registration/reference_files/mask.mnc \
	--ncontrols 10 \
	--threshold 5 \
	--dataset $dataset \
	--parallel true \
	--nproc 4
    
    
# Downsample images ----------------------------------------------------------
    
echo "Resampling template and mask images..."

#Downsample model file to 3.0mm
autocrop -quiet -clobber -isostep 3.0 \
	data/human/registration/reference_files/model.mnc \
	data/human/registration/reference_files/model_3.0mm.mnc

#Downsample mask file to 3.0mm
autocrop -quiet -clobber -isostep 3.0 \
	data/human/registration/reference_files/mask.mnc \
	data/human/registration/reference_files/mask_3.0mm.mnc
    
#Downsample model file to 1.0mm
autocrop -quiet -clobber -isostep 1.0 \
	data/human/registration/reference_files/model.mnc \
	data/human/registration/reference_files/model_1.0mm.mnc

#Downsample mask file to 1.0mm
autocrop -quiet -clobber -isostep 1.0 \
	data/human/registration/reference_files/mask.mnc \
	data/human/registration/reference_files/mask_1.0mm.mnc
    
echo "Resampling absolute effect size images..."

#Downsample absolute effect size images to 3.0mm
python3 resample_minc_images.py \
	--imgdir ${human_dir}/effect_sizes/absolute/resolution_0.5/ \
	--outdir ${human_dir}/effect_sizes/absolute/resolution_3.0/ \
	--isostep 3.0 \
	--parallel true \
	--nproc 4

#Downsample absolute effect size images to 1.0mm
python3 resample_minc_images.py \
	--imgdir ${human_dir}/effect_sizes/absolute/resolution_0.5/ \
	--outdir ${human_dir}/effect_sizes/absolute/resolution_1.0/ \
	--isostep 1.0 \
	--parallel true \
	--nproc 4
    
echo "Resampling relative effect size images..."

#Downsample relative effect size images to 3.0mm
python3 resample_minc_images.py \
	--imgdir ${human_dir}/effect_sizes/relative/resolution_0.5/ \
	--outdir ${human_dir}/effect_sizes/relative/resolution_3.0/ \
	--isostep 3.0 \
	--parallel true \
	--nproc 4

#Downsample relative effect size images to 1.0mm
python3 resample_minc_images.py \
	--imgdir ${human_dir}/effect_sizes/relative/resolution_0.5/ \
	--outdir ${human_dir}/effect_sizes/relative/resolution_1.0/ \
	--isostep 1.0 \
	--parallel true \
	--nproc 4
    
    
# Construct effect sizes matrices --------------------------------------------
    
echo "Build absolute effect size matrix..."

#Create matrix using absolute effect sizes at 3.0mm
python3 build_effect_size_matrix.py \
	--imgdir ${human_dir}/effect_sizes/absolute/resolution_3.0/ \
	--outfile ${human_dir}/effect_sizes/absolute/human_effect_sizes_absolute_3.0mm.csv \
	--mask data/human/registration/reference_files/mask_3.0mm.mnc
    
#Create matrix using absolute effect sizes at 1.0mm
python3 build_effect_size_matrix.py \
	--imgdir ${human_dir}/effect_sizes/absolute/resolution_1.0/ \
	--outfile ${human_dir}/effect_sizes/absolute/human_effect_sizes_absolute_1.0mm.csv \
	--mask data/human/registration/reference_files/mask_1.0mm.mnc 
    
echo "Build relative effect size matrix..."    

#Create matrix using absolute effect sizes at 3.0mm
python3 build_effect_size_matrix.py \
	--imgdir ${human_dir}/effect_sizes/relative/resolution_3.0/ \
	--outfile ${human_dir}/effect_sizes/relative/human_effect_sizes_relative_3.0mm.csv \
	--mask data/human/registration/reference_files/mask_3.0mm.mnc

#Create matrix using absolute effect sizes at 1.0mm
python3 build_effect_size_matrix.py \
	--imgdir ${human_dir}/effect_sizes/relative/resolution_1.0/ \
	--outfile ${human_dir}/effect_sizes/relative/human_effect_sizes_relative_1.0mm.csv \
	--mask data/human/registration/reference_files/mask_1.0mm.mnc
    
    
# Run clustering -------------------------------------------------------------
    
echo "Clustering human effect size images..."
    
#Cluster human effect sizes
Rscript cluster_human_data.R \
	--file1 ${human_dir}/effect_sizes/absolute/human_effect_sizes_absolute_3.0mm.csv \
	--file2 ${human_dir}/effect_sizes/relative/human_effect_sizes_relative_3.0mm.csv \
	--rownames file \
	--nclusters 10 \
	--K 10 \
	--t 20 \
	--sigma 0.5 \
	--metric cor \
	--outfile ${human_dir}/clustering/human_clusters_groups10_3.0mm.csv


# Generate representative cluster maps ---------------------------------------

echo "Creating cluster maps using absolute effect sizes and mean aggregation..."

# Create cluster maps using absolute effect sizes at 3.0mm
Rscript create_cluster_maps.R \
	--clusterfile ${human_dir}/clustering/human_clusters_groups10_3.0mm.csv \
	--imgdir ${human_dir}/effect_sizes/absolute/resolution_3.0/ \
	--outdir ${human_dir}/clustering/cluster_maps/absolute/resolution_3.0/mean/ \
	--method mean \
	--jacobians absolute
    
# Create cluster maps using absolute effect sizes at 1.0mm
Rscript create_cluster_maps.R \
	--clusterfile ${human_dir}/clustering/human_clusters_groups10_3.0mm.csv \
	--imgdir ${human_dir}/effect_sizes/absolute/resolution_1.0/ \
	--outdir ${human_dir}/clustering/cluster_maps/absolute/resolution_1.0/mean/ \
	--method mean \
	--jacobians absolute

echo "Creating cluster maps using relative effect sizes and mean aggregation..."

# Create cluster maps using relative effect sizes at 3.0mm
Rscript create_cluster_maps.R \
	--clusterfile ${human_dir}/clustering/human_clusters_groups10_3.0mm.csv \
	--imgdir ${human_dir}/effect_sizes/relative/resolution_3.0/ \
	--outdir ${human_dir}/clustering/cluster_maps/relative/resolution_3.0/mean/ \
	--method mean \
	--jacobians relative
    
# Create cluster maps using relative effect sizes at 1.0mm
Rscript create_cluster_maps.R \
	--clusterfile ${human_dir}/clustering/human_clusters_groups10_3.0mm.csv \
	--imgdir ${human_dir}/effect_sizes/relative/resolution_1.0/ \
	--outdir ${human_dir}/clustering/cluster_maps/relative/resolution_1.0/mean/ \
	--method mean \
	--jacobians relative
    
    
deactivate
