#!/bin/bash

module purge
module load \
    zlib \
	minc-toolkit/1.9.18.1-mp7vcse \
	r/3.6.3-zlk4uk6 \
	r-packages/2022-01-26 \
	ants

source ../../.venv/bin/activate

#Download and save world coordinates as .csv file
Rscript ../../fetch_microarray_coordinates.R \
    --metadata data/SampleInformation_pipeline_v1.csv \
    --outfile data/AHBA_microarray_coordinates_mni.csv \
    --labels true

#Create a mask image for voxels corresponding to coordinates
Rscript ../../coordinates_to_minc.R \
     --coordinates data/AHBA_microarray_coordinates_mni.csv \
     --template data/mni_icbm152_t1_tal_nlin_sym_09c_extracted.mnc \
     --outfile data/AHBA_microarray_mask_mni.mnc \
     --type mask
     
Rscript ../../coordinates_to_minc.R \
     --coordinates data/AHBA_microarray_coordinates_mni.csv \
     --template data/mni_icbm152_t1_tal_nlin_sym_09c_extracted.mnc \
     --outfile data/AHBA_microarray_labels_mni_tmp.mnc \
     --type labels
     
deactivate

module load minc-stuffs

minc_label_ops --convert \
    data/AHBA_microarray_labels_mni_tmp.mnc \
    data/AHBA_microarray_labels_mni.mnc
rm data/AHBA_microarray_labels_mni_tmp.mnc

module purge

#Examine alignment
#register \
#    mni_icbm152_t1_tal_nlin_sym_09c_extracted.mnc \
#    AHBA_microarray_mask_mni.mnc
    
# -----------------------
# AHBA_microarray_mask_mni.mnc is correcty aligned with mni_icbm152_t1_tal_nlin_sym_09c_extracted.mnc
# What if we transform the mask to the consensus average space? 
# -----------------------

module load \
    zlib \
	minc-toolkit/1.9.18.1-mp7vcse \
	r/3.6.3-zlk4uk6 \
	r-packages/2022-01-26 \
	ants

source ../../.venv/bin/activate

#Transform the coordinates mask from MNI space to study space 
antsApplyTransforms \
	-d 3 \
	--input data/AHBA_microarray_mask_mni.mnc \
	--output data/AHBA_microarray_mask_mni_to_study_tmp.mnc \
	--reference-image data/template_sharpen_shapeupdate_1.0mm.mnc \
	-t [transforms/average_to_icbm_nlin_sym_09c_minc0_GenericAffine.xfm,1] \
	-t transforms/average_to_icbm_nlin_sym_09c_minc1_inverse_NL.xfm
    
deactivate 

module load minc-stuffs

minc_label_ops --convert \
    data/AHBA_microarray_mask_mni_to_study_tmp.mnc \
    data/AHBA_microarray_mask_mni_to_study.mnc
rm data/AHBA_microarray_mask_mni_to_study_tmp.mnc

module purge

#Examine alignment
#register \
#    template_sharpen_shapeupdate_1.0mm.mnc \
#    AHBA_microarray_mask_mni_to_study.mnc

# -----------------------
# AHBA_microarray_mask_mni_to_study.mnc is correcty aligned with template_sharpen_shapeupdate.mnc
#
# So the masks work, but we lose sample information in doing this. What happens if we naively try 
# to transform the coordinates themselves?
# -----------------------

module load \
    zlib \
	minc-toolkit/1.9.18.1-mp7vcse \
	r/3.6.3-zlk4uk6 \
	r-packages/2022-01-26 \
	ants

source ../../.venv/bin/activate

antsApplyTransformsToPoints \
   -d 3 \
   --input data/AHBA_microarray_coordinates_mni.csv \
   --output data/AHBA_microarray_coordinates_mni_to_study.csv \
	-t [transforms/average_to_icbm_nlin_sym_09c_minc0_GenericAffine.xfm,1] \
	-t transforms/average_to_icbm_nlin_sym_09c_minc1_inverse_NL.xfm \
   --precision 1
    
Rscript ../../coordinates_to_minc.R \
   --coordinates data/AHBA_microarray_coordinates_mni_to_study.csv \
   --template data/template_sharpen_shapeupdate_1.0mm.mnc \
   --outfile data/AHBA_microarray_mask_study.mnc \
   --type mask
   
#Examine alignment
#register \
#    template_sharpen_shapeupdate_1.0mm.mnc \
#    AHBA_microarray_mask_study.mnc
    
# ------------------
# AHBA_microarray_mask_study.mnc is not aligned with template_sharpen_shapeupdate.mnc
#
# Ben said that coordinates use the inverse transforms from images. So let's try that. 
# -------------------

antsApplyTransformsToPoints \
   -d 3 \
   --input data/AHBA_microarray_coordinates_mni.csv \
   --output data/AHBA_microarray_coordinates_mni_to_study_using_inverses.csv \
	-t [transforms/average_to_icbm_nlin_sym_09c_minc0_GenericAffine.xfm,0] \
	-t transforms/average_to_icbm_nlin_sym_09c_minc1_NL.xfm \
   --precision 1
    
Rscript ../../coordinates_to_minc.R \
   --coordinates data/AHBA_microarray_coordinates_mni_to_study_using_inverses.csv \
   --template data/template_sharpen_shapeupdate_1.0mm.mnc \
   --outfile data/AHBA_microarray_mask_study_using_inverses.mnc \
   --type mask
   
#Examine alignment
#register \
#    template_sharpen_shapeupdate_1.0mm.mnc \
#    AHBA_microarray_mask_study_using_inverses.mnc

# ------------------
# It looks like this might actually have worked. The samples are at least within the brain area.
#
# Let's create a pseudo-atlas for these samples
# -------------------

Rscript ../../coordinates_to_minc.R \
   --coordinates data/AHBA_microarray_coordinates_mni_to_study_using_inverses.csv \
   --template data/template_sharpen_shapeupdate_1.0mm.mnc \
   --outfile data/AHBA_microarray_labels_study_using_inverses_tmp.mnc \
   --type labels

deactivate

module load minc-stuffs

minc_label_ops --convert \
    data/AHBA_microarray_labels_study_using_inverses_tmp.mnc \
    data/AHBA_microarray_labels_study_using_inverses.mnc
rm data/AHBA_microarray_labels_study_using_inverses_tmp.mnc

module purge
