#!/bin/bash

Rscript normative_growth_normalization.R \
    --demographics data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc.csv \
    --voxels data/human/derivatives/POND_SickKids/combat/jacobians_normalized.csv \
    --outfile data/human/derivatives/POND_SickKids/effect_sizes/absolute/tmp/tmp.csv \
    --key file \
    --df 5 \
    --combat true \
    --combat-batch Site \
    --parallel true \
    --nproc 8
    
    
# Rscript calculate_human_effect_sizes.R \
#     --demographics data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc.csv \
#     --maskfile data/human/registration/reference_files/mask.mnc \
#     --imgdir data/human/derivatives/POND_SickKids/jacobians/absolute \
#     --outdir tmp/ \
#     --ncontrols 10 \
#     --parallel true \
#     --nproc 4