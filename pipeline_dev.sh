#!/bin/bash

source activate_venv.sh

# Get mouse clustering files from Jacob's directories
# source get_mouse_files.sh


# Generate mouse and human cluster masks -------------------------------------

# Create cluster masks using multiple thresholds
thresholds=(0.1 0.2 0.3 0.4 0.5 0.6 0.8 0.7 0.9 1.0)
for i in ${thresholds[@]};
do
    echo "Creating masks at threshold ${i} ..."

    # Create mouse mask using absolute effect sizes at 200um
    python3 create_cluster_masks.py \
	    --imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
	    --outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_$i/ \
	    --threshold $i \
	    --symmetric true
        
    # Create mouse mask using relative effect sizes at 200um    
    python3 create_cluster_masks.py \
	    --imgdir data/mouse/clustering/cluster_maps/relative/resolution_200/mean/ \
	    --outdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_$i/ \
	    --threshold $i \
	    --symmetric true

    # Create human mask using absolute effect sizes at 3.0mm
    python3 create_cluster_masks.py \
	    --imgdir data/human/clustering/cluster_maps/absolute/resolution_3.0/mean/ \
	    --outdir data/human/clustering/cluster_masks/absolute/resolution_3.0/mean/threshold_$i/ \
	    --threshold $i \
	    --symmetric true

    # Create human mask using relative effect sizes at 3.0mm        
    python3 create_cluster_masks.py \
	    --imgdir data/human/clustering/cluster_maps/relative/resolution_3.0/mean/ \
	    --outdir data/human/clustering/cluster_masks/relative/resolution_3.0/mean/threshold_$i/ \
	    --threshold $i \
	    --symmetric true
        
    # Create human mask using absolute effect sizes at 1.0mm
    python3 create_cluster_masks.py \
	    --imgdir data/human/clustering/cluster_maps/absolute/resolution_1.0/mean/ \
	    --outdir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_$i/ \
	    --threshold $i \
	    --symmetric true

    # Create human mask using relative effect sizes at 1.0mm
    python3 create_cluster_masks.py \
	    --imgdir data/human/clustering/cluster_maps/relative/resolution_1.0/mean/ \
	    --outdir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_$i/ \
	    --threshold $i \
	    --symmetric true

done


# Re-organize gene expression files ------------------------------------------

#echo "Organizing gene expression files..."
#source organize_expression_files.sh


# Transform human microarray samples -----------------------------------------

echo "Mapping AHBA microarray samples to imaging study space..."

#Download microarray sample coordinates in MNI space
Rscript fetch_microarray_coordinates.R \
	--metadata data/human/expression/SampleInformation_pipeline_abagen.csv \
	--outfile data/human/expression/AHBA_microarray_coordinates_mni.csv \
	--labels true
    
#Transform microarray sample coordinates from MNI space to study space
#Coordinates need to be transformed using inverse transforms
antsApplyTransformsToPoints \
	-d 3 \
	--input data/human/expression/AHBA_microarray_coordinates_mni.csv \
	--output data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
	-t [data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc0_GenericAffine.xfm,0] \
	-t data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc1_NL.xfm
   
# Create a mask indicating sample positions
Rscript coordinates_to_minc.R \
	--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
	--template data/human/registration/reference_files/model_1.0mm.mnc \
	--outfile data/human/expression/AHBA_microarray_mask_studyspace_1.0mm.mnc \
	--type mask
    
# Create a pseudo-atlas for individual samples
Rscript coordinates_to_minc.R \
	--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
	--template data/human/registration/reference_files/model_1.0mm.mnc \
	--outfile data/human/expression/AHBA_microarray_labels_studyspace_1.0mm.mnc \
	--type labels
    

# Generate mouse and human cluster signatures using input gene space ---------

# Generate cluster signatures using multiple thresholds
thresholds=(0.1 0.2 0.3 0.4 0.5 0.6 0.8 0.7 0.9 1.0)
for i in ${thresholds[@]};
do
	echo "Creating mouse cluster signatures using threshold $i masks and input gene space..."
    
    # Mouse cluster signatures using absolute masks at 200um
    python3 mouse_cluster_signatures.py \
	    --clusterdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_$i/ \
	    --exprdir data/mouse/expression/input_space/ \
	    --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	    --outfile data/mouse/cluster_signatures/input_space/mouse_cluster_signatures_abs_mean_threshold${i}_inputspace.csv \
	    --parallel true \
	    --nproc 4

    # Mouse cluster signatures using relative masks at 200um        
    python3 mouse_cluster_signatures.py \
	    --clusterdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_$i/ \
	    --exprdir data/mouse/expression/input_space/ \
	    --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	    --outfile data/mouse/cluster_signatures/input_space/mouse_cluster_signatures_rel_mean_threshold${i}_inputspace.csv \
	    --parallel true \
	    --nproc 4
        
    echo "Creating human cluster signatures using threshold $i masks and input gene space..."

    # Human cluster signatures using absolute masks at 1.0mm
    Rscript human_cluster_signatures.R \
	    --clusterdir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_$i/ \
	    --exprdir data/human/expression/input_space/ \
	    --coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
	    --template data/human/registration/reference_files/model_1.0mm.mnc \
	    --outfile data/human/cluster_signatures/input_space/human_cluster_signatures_abs_mean_threshold${i}_inputspace.csv \
	    --parallel true \
	    --nproc 4
    
    # Human cluster signatures using relative masks at 1.0mm    
    Rscript human_cluster_signatures.R \
	    --clusterdir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_$i/ \
	    --exprdir data/human/expression/input_space/ \
	    --coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
	    --template data/human/registration/reference_files/model_1.0mm.mnc \
	    --outfile data/human/cluster_signatures/input_space/human_cluster_signatures_rel_mean_threshold${i}_inputspace.csv \
	    --parallel true \
	    --nproc 4
        
done


# Generate mouse and human cluster signatures using latent spaces ------------

# for i in ${thresholds[@]};
# do
        
#     echo "Creating mouse cluster signatures using threshold $i masks and latent gene space..."

#     python3 mouse_cluster_signatures.py \
#         --clusterdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_$i/ \
#         --exprdir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold${i}_latentspace100.csv \
#         --parallel true \
#         --nproc 8
        
#     python3 mouse_cluster_signatures.py \
#         --clusterdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_$i/ \
#         --exprdir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold${i}_latentspace100.csv \
#         --parallel true \
#         --nproc 8
    
#     echo "Creating human cluster signatures using threshold $i masks and latent gene space..."

#     Rscript human_cluster_signatures.R \
#         --clusterdir data/human/clustering/cluster_masks/absolute/resolution_3.0/mean/threshold_$i/ \
#         --exprdir data/human/expression/latent_space_100/ \
#         --metadata data/human/SampleInformation_pipeline_v1.csv \
#         --template data/human/registration/reference_files/model_3.0mm.mnc \
#         --outfile data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold${i}_latentspace100.csv \
#         --parallel true \
#         --nproc 8
        
#     Rscript human_cluster_signatures.R \
#         --clusterdir data/human/clustering/cluster_masks/relative/resolution_3.0/mean/threshold_$i/ \
#         --exprdir data/human/expression/latent_space_100/ \
#         --metadata data/human/SampleInformation_pipeline_v1.csv \
#         --template data/human/registration/reference_files/model_3.0mm.mnc \
#         --outfile data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold${i}_latentspace100.csv \
#         --parallel true \
#         --nproc 8
        
# done




# for i in ${thresholds[@]};
# do

#     echo "Computing latent space similarity matrices using threshold $i data..."

#     echo "Using absolute jacobians..."
#     Rscript latent_space_similarity_matrix.R \
#         --mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold${i}_latentspace100.csv \
#         --human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold${i}_latentspace100.csv \
#         --metric correlation \
#         --outfile data/similarity/latent_space_100/similarity_hm_abs_mean_threshold${i}_latentspace100.csv \
#         --save-intermediate false

#     echo "Using relative jacobians..."
#     Rscript latent_space_similarity_matrix.R \
#         --mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold${i}_latentspace100.csv \
#         --human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold${i}_latentspace100.csv \
#         --metric correlation \
#         --outfile data/similarity/latent_space_100/similarity_hm_rel_mean_threshold${i}_latentspace100.csv \
#         --save-intermediate false

# done

deactivate
