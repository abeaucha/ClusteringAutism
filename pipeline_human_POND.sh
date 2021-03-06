#!/bin/bash

source activate_venv

#Extract compressed jacobians and convert to minc
#python3 extract_zipped_images.py \
#	--imgdir data/human/registration/jacobians/absolute/smooth/ \
#	--parallel true \
#	--nproc 8

#python3 extract_zipped_images.py \
#	--imgdir data/human/registration/jacobians/relative/smooth/ \
#	--parallel true \
#	--nproc 8

echo "Calculating human effect sizes using absolute jacobians..."

#Calculate human effect sizes based on relative and absolute volumes
Rscript calculate_human_effect_sizes.R \
	--demographics data/human/registration/DBM_input_demo_passedqc.csv \
	--imgdir data/human/registration/jacobians/absolute/smooth/ \
	--outdir data/human/POND/effect_sizes/absolute/resolution_0.5/ \
	--maskfile data/human/registration/reference_files/mask.mnc \
	--ncontrols 10 \
	--dataset POND \
	--parallel true \
	--nproc 4
    
echo "Calculating human effect sizes using relative jacobians..."

Rscript calculate_human_effect_sizes.R \
	--demographics data/human/registration/DBM_input_demo_passedqc.csv \
	--imgdir data/human/registration/jacobians/relative/smooth/ \
	--outdir data/human/POND/effect_sizes/relative/resolution_0.5/ \
	--maskfile data/human/registration/reference_files/mask.mnc \
	--ncontrols 10 \
	--dataset POND \
	--parallel true \
	--nproc 4
    
echo "Resampling template and mask images..."

#Downsample model and mask files
autocrop -quiet -clobber -isostep 3.0 \
	data/human/registration/reference_files/model.mnc \
	data/human/registration/reference_files/model_3.0mm.mnc

autocrop -quiet -clobber -isostep 3.0 \
	data/human/registration/reference_files/mask.mnc \
	data/human/registration/reference_files/mask_3.0mm.mnc

echo "Resampling absolute effect size images..."

#Downsample human effect size images 
python3 resample_minc_images.py \
	--imgdir data/human/POND/effect_sizes/absolute/resolution_0.5/ \
	--outdir data/human/POND/effect_sizes/absolute/resolution_3.0/ \
	--isostep 3.0 \
	--parallel true \
	--nproc 4
    
echo "Resampling relative effect size images..."

python3 resample_minc_images.py \
	--imgdir data/human/POND/effect_sizes/relative/resolution_0.5/ \
	--outdir data/human/POND/effect_sizes/relative/resolution_3.0/ \
	--isostep 3.0 \
	--parallel true \
	--nproc 4
    
echo "Build absolute effect size matrix..."

#Create matrices for human effect sizes
python3 build_effect_size_matrix.py \
	--imgdir data/human/POND/effect_sizes/absolute/resolution_3.0/ \
	--outfile data/human/POND/effect_sizes/absolute/POND_effect_sizes_absolute_3.0mm.csv \
	--mask data/human/registration/reference_files/mask_3.0mm.mnc
    
echo "Build relative effect size matrix..."    

python3 build_effect_size_matrix.py \
	--imgdir data/human/POND/effect_sizes/relative/resolution_3.0/ \
	--outfile data/human/POND/effect_sizes/relative/POND_effect_sizes_relative_3.0mm.csv \
	--mask data/human/registration/reference_files/mask_3.0mm.mnc
    
echo "Clustering human effect size images..."
    
#Cluster human effect sizes
Rscript cluster_human_data.R \
	--file1 data/human/POND/effect_sizes/absolute/POND_effect_sizes_absolute_3.0mm.csv \
	--file2 data/human/POND/effect_sizes/relative/POND_effect_sizes_relative_3.0mm.csv \
	--rownames file \
	--nclusters 10 \
	--K 10 \
	--t 20 \
	--sigma 0.5 \
	--metric cor \
	--outfile data/human/POND/clustering/POND_clusters_groups10_3.0mm.csv
    
echo "Creating cluster maps using absolute effect sizes and mean aggregation..."

# Create cluster maps for abs and rel volumes using mean and median
Rscript create_cluster_maps.R \
	--clusterfile data/human/POND/clustering/POND_clusters_groups10_3.0mm.csv \
	--imgdir data/human/POND/effect_sizes/absolute/resolution_3.0/ \
	--outdir data/human/POND/clustering/cluster_maps/absolute/resolution_3.0/mean/ \
	--method mean \
	--jacobians absolute
    
echo "Creating cluster maps using relative effect sizes and mean aggregation..."

Rscript create_cluster_maps.R \
	--clusterfile data/human/POND/clustering/POND_clusters_groups10_3.0mm.csv \
	--imgdir data/human/POND/effect_sizes/relative/resolution_3.0/ \
	--outdir data/human/POND/clustering/cluster_maps/relative/resolution_3.0/mean/ \
	--method mean \
	--jacobians relative
    
deactivate
