#!/bin/bash -l

source activate_venv.sh

#Input files
training=data/mouse/expression/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_GM_labelled_67_scaled.csv
mouse_input=data/mouse/expression/input_space/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv
human_input=data/human/expression/input_space/HumanExpressionMatrix_samples_pipeline_abagen_homologs_scaled.csv

#Output directories
outdir=data/MLP_outcomes/
if [ ! -d "$outdir" ]; then
	mkdir -p "$outdir"
fi

mouse_latent_space=data/mouse/expression/latent_space/
if [ ! -d "$mouse_latent_space" ]; then
	mkdir -p "$mouse_latent_space"
fi

human_latent_space=data/human/expression/latent_space/
if [ ! -d "$human_latent_space" ]; then
	mkdir -p "$human_latent_space"
fi

#Parameters file
current_date=$(date +'%Y%m%d')
paramsfile=${outdir}MLP_params_${current_date}.csv
touch $paramsfile

niterations=500
nunits=200
L2=0.0
nepochs=150
totalsteps=200
learningrate=1e-05
optimizer=AdamW
confusionmatrix=false

for i in $(seq 1 $niterations);
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
		--training $training \
		--mousetransform $mouse_input \
		--humantransform $human_input \
		--nunits $nunits \
		--L2 $L2 \
		--nepochs $nepochs \
		--totalsteps $totalsteps \
		--learningrate $learningrate \
        --optimizer $optimizer \
		--confusionmatrix false \
		--seed $i \
		--saveparams true \
		--paramsheader $paramsheader \
		--verbose true

	cat ${outdir}MLP_params.csv >> $paramsfile

	mv ${outdir}MLP_labels67_layers3_units${nunits}_L2${L2}_humanprob.csv \
		${outdir}MLP_labels67_layers3_units${nunits}_L2${L2}_humanprob_$i.csv

	mv ${outdir}MLP_labels67_layers3_units${nunits}_L2${L2}_humantransform.csv \
		${human_latent_space}MLP_labels67_layers3_units${nunits}_L2${L2}_humantransform_$i.csv

	mv ${outdir}MLP_labels67_layers3_units${nunits}_L2${L2}_mouseprob.csv \
		${outdir}MLP_labels67_layers3_units${nunits}_L2${L2}_mouseprob_$i.csv

	mv ${outdir}MLP_labels67_layers3_units${nunits}_L2${L2}_mousetransform.csv \
		${mouse_latent_space}MLP_labels67_layers3_units${nunits}_L2${L2}_mousetransform_$i.csv

done

rm ${outdir}MLP_params.csv

deactivate
