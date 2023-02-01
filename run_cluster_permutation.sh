#!/bin/bash

# run_cluster_permutation.sh
# Author: Antoine Beauchamp


#Activate virtual environment
source activate_venv.sh

#Define some variables
outdir=data/human/permutation/
nk=5
np=200
jacobians=(absolute relative)
jacobians_short=(abs rel)
signs=(positive negative)

echo "Permutating cluster labels..."
Rscript permute_cluster_labels.R \
    --infile data/human/clustering/human_clusters_groups10_3.0mm.csv \
    --outdir ${outdir}clusters/ \
    --nk $nk \
    --np $np

echo "Running permutations..."
for p in $(seq 1 $np);
do

    echo "On permutation $p of $np"

    clusterfile=${outdir}clusters/human_clusters_groups10_3.0mm_nk_5_permutation_${p}.csv
    outdir_p=${outdir}/permutation_${p}/

    #Iterate over jacobians
    for j in $(seq 1 ${#jacobians[@]});
    do
        
        jacobian=${jacobians[$j-1]}
        jacobian_short=${jacobians_short[$j-1]}
    
        #Create cluster maps using effect sizes at 1.0mm
        echo "Creating $jacobian cluster maps..."
        Rscript create_cluster_maps.R \
            --clusterfile $clusterfile \
            --imgdir data/human/effect_sizes/${jacobian}/resolution_1.0/ \
            --outdir ${outdir_p}/cluster_maps/${jacobian}/ \
            --method mean \
            --jacobians ${jacobian}
            
        #Create cluster masks with thresholding at 0.2
        echo "Creating ${jacobian} cluster masks..."
        python3 create_image_masks.py \
            --imgdir ${outdir_p}/cluster_maps/${jacobian}/ \
            --outdir ${outdir_p}/cluster_masks/${jacobian}/ \
            --mask data/human/registration/reference_files/mask_1.0mm.mnc \
            --method top_n \
            --threshold 0.2 \
            --symmetric true \
            --signed true
            
        #Create cluster signatures and similarity matrices
        #Iterate over mask signs
        for s in $(seq 1 ${#signs[@]});
        do
        
            sign=${signs[$s-1]}
        
            echo "Creating ${jacobian} cluster signatures using ${sign} mask..."
            Rscript human_cluster_signatures.R \
                --cluster-dir ${outdir_p}/cluster_masks/${jacobian}/ \
                --expr-dir data/human/expression/latent_space_100/ \
                --coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
                --template data/human/registration/reference_files/model_1.0mm.mnc \
                --sign ${sign} \
                --outfile ${outdir_p}/cluster_signatures/human_cluster_signatures_${jacobian_short}_${sign}.csv \
                --parallel true \
                --nproc 8
                    
            echo "Creating similarity matrix using ${jacobian} ${sign} mask signatures..."
            Rscript latent_space_similarity_matrix.R \
                --mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_${jacobian_short}_mean_threshold_topn_0.2_${sign}_latentspace100.csv \
                --human ${outdir_p}/cluster_signatures/human_cluster_signatures_${jacobian_short}_${sign}.csv \
                --metric correlation \
                --outfile ${outdir_p}/similarity_matrix/similarity_hm_${jacobian_short}_${sign}.csv \
                --save-intermediate false
                
        done
    done
done
    
deactivate