#!/bin/bash


# Calculate effect sizes -----------------------------------------------------

#Calculate human effect sizes based using absolute volumes
echo "Calculating human effect sizes using absolute jacobians..."
Rscript calculate_human_effect_sizes.R \
	--demographics data/human/registration/DBM_input_demo_passedqc.csv \
	--imgdir data/human/registration/jacobians/absolute/smooth/ \
	--outdir data/human/effect_sizes/absolute/resolution_0.5/ \
	--maskfile data/human/registration/reference_files/mask.mnc \
	--ncontrols 10 \
	--threshold 5 \
	--parallel true \
	--nproc 4