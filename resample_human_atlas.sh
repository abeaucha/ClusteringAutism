#!/bin/bash

module load minc-stuffs

antsApplyTransforms \
	-d 3 \
	--input data/human/atlas/annotation_full.mnc \
	--output data/human/atlas/annotation_full_studyspace_tmp.mnc \
	--reference-image data/human/registration/average_to_MNI/template_sharpen_shapeupdate_1.0mm.mnc \
	--interpolation NearestNeighbor \
	-t [data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09b_minc0_GenericAffine.xfm,1] \
	-t data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09b_minc1_inverse_NL.xfm 
 
minc_label_ops --convert \
data/human/atlas/annotation_full_studyspace_tmp.mnc \
data/human/atlas/annotation_full_studyspace.mnc



# antsApplyTransformsToPoints \
# 	-d 3 \
# 	--input data/human/expression/AHBA_microarray_coordinates_mni.csv \
# 	--output data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
# 	-t [data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc0_GenericAffine.xfm,0] \
# 	-t data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc1_NL.xfm

