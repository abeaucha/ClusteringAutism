#!/bin/bash

outdir=data/mouse/expression/
if [ ! -d '$outdir' ];
then
    mkdir -p $outdir
fi

mouse_expr_dir=/projects/abeauchamp/Projects/MouseHumanMapping/Paper_TranscriptomicSimilarity/main/AMBA/data/

cp ${mouse_expr_dir}MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed.csv ${outdir}MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed.csv

cp ${mouse_expr_dir}MouseExpressionTree_DSURQE.RData ${outdir}


outdir=data/human/expression/
if [ ! -d '$outdir' ];
then
    mkdir -p $outdir
fi

human_expr_dir=/projects/abeauchamp/Projects/MouseHumanMapping/Paper_TranscriptomicSimilarity/main/AHBA/data/

#Filenames subject to change
cp ${human_expr_dir}HumanExpressionMatrix_samples_pipeline_v1.csv ${outdir}
cp ${human_expr_dir}SampleInformation_pipeline_v1.csv ${outdir}
cp ${human_expr_dir}HumanExpressionTree.RData ${outdir}

expr_dir=/projects/abeauchamp/Projects/MouseHumanMapping/Paper_TranscriptomicSimilarity/main/data/

cp ${expr_dir}MouseHumanGeneHomologs.csv data/
cp ${expr_dir}TreeLabels.RData data/
cp ${expr_dir}TreeLabelsReordered.RData data/
