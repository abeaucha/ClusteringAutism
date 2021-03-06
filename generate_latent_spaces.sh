#!/bin/bash -l

source activate_venv

outdir=data/MLP_outcomes/

if [ ! -d "$outdir" ]; then
	mkdir "$outdir"
fi

current_date=$(date +'%Y%m%d')
paramsfile=${outdir}MLP_params_${current_date}.csv
touch $paramsfile

for i in {1..500}
do
 echo "Iteration $i"

 if [ $i -eq 1 ];
 then
	 paramsheader=true
 else
	 paramsheader=false
 fi

 python3 train_multilayer_perceptron.py \
	 --outdir $outdir \
	 --training data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_GM_labelled_67_scaled.csv \
	 --mousetransform data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv \
	 --humantransform data/HumanExpressionMatrix_samples_pipeline_v1_homologs_scaled.csv \
	 --nunits 200 \
	 --L2 1e-6 \
	 --nepochs 200 \
	 --learningrate 1e-5 \
	 --confusionmatrix true \
	 --seed $i \
	 --saveparams true \
	 --paramsheader $paramsheader \
	 --verbose true

 cat ${outdir}MLP_params.csv >> $paramsfile

 mv ${outdir}MLP_labels67_layers3_units200_L21e-06_confusionmatrix_training.csv \
	 ${outdir}MLP_labels67_layers3_units200_L21e-06_confusionmatrix_training_$i.csv

 mv ${outdir}MLP_labels67_layers3_units200_L21e-06_humanprob.csv \
	 ${outdir}MLP_labels67_layers3_units200_L21e-06_humanprob_$i.csv

 mv ${outdir}MLP_labels67_layers3_units200_L21e-06_humantransform.csv \
	 ${outdir}MLP_labels67_layers3_units200_L21e-06_humantransform_$i.csv

 mv ${outdir}MLP_labels67_layers3_units200_L21e-06_mouseprob.csv \
	 ${outdir}MLP_labels67_layers3_units200_L21e-06_mouseprob_$i.csv

 mv ${outdir}MLP_labels67_layers3_units200_L21e-06_mousetransform.csv \
	 ${outdir}MLP_labels67_layers3_units200_L21e-06_mousetransform_$i.csv

done

rm ${outdir}MLP_params.csv

deactivate
