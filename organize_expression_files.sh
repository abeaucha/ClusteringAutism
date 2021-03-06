#!/bin/bash

#Create mouse expression directory
if [ ! -d "data/mouse/expression/" ]; then
    mkdir data/mouse/expression
fi

#Create mouse input space directory
if [ ! -d "data/mouse/expression/input_space/" ]; then
    mkdir data/mouse/expression/input_space/
fi

#Create symlink for input space data
cd data/mouse/expression/input_space/
input_space=MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv
if [ ! -f $input_space ];
then
    ln -s ../../../$input_space ./
fi

#Return to main directory
cd ../../../..

#Create mouse latent space directory
if [ ! -d "data/mouse/expression/latent_space/" ]; then
    mkdir data/mouse/expression/latent_space/
fi

#Create symlinks for mouse latent space data
cd data/mouse/expression/latent_space/
mouse_latent_spaces=(../../../MLP_outcomes/*mousetransform*.csv)
for file in ${mouse_latent_spaces[@]}; 
do
    filename=$(basename -- "$file")
    if [ ! -f $filename ];
    then
        ln -s $file ./
    fi
done

#Return to main directory
cd ../../../..

if [ ! -d 'data/mouse/expression/latent_space_100/' ];
then
    mkdir data/mouse/expression/latent_space_100
fi

#Create symlinks for mouse latent space data
cd data/mouse/expression/latent_space_100/
mouse_latent_spaces=(../../../MLP_outcomes/*mousetransform_1[0-9][0-9].csv)
for file in ${mouse_latent_spaces[@]}; 
do
    filename=$(basename -- "$file")
    if [ ! -f $filename ];
    then
        ln -s $file ./
    fi
done

#Return to main directory
cd ../../../..

#Create human expression directory
if [ ! -d "data/human/expression/" ]; then
    mkdir data/human/expression
fi

#Create human input space directory
if [ ! -d "data/human/expression/input_space/" ]; then
    mkdir data/human/expression/input_space/
fi

# #Create symlink for input space data
cd data/human/expression/input_space/
input_space=HumanExpressionMatrix_samples_pipeline_v1_homologs_scaled.csv
if [ ! -f $input_space ];
then
    ln -s ../../../$input_space ./
fi

#Return to main directory
cd ../../../..

#Create human latent space directory
if [ ! -d "data/human/expression/latent_space/" ]; then
    mkdir data/human/expression/latent_space/
fi

#Create symlinks for human latent space data
cd data/human/expression/latent_space/
human_latent_spaces=(../../../MLP_outcomes/*humantransform*.csv)
for file in ${human_latent_spaces[@]}; 
do
    filename=$(basename -- "$file")
    if [ ! -f $filename ];
    then
        ln -s $file ./
    fi
done

#Return to main directory
cd ../../../..

#Create human latent space directory
if [ ! -d "data/human/expression/latent_space_100/" ]; 
then
    mkdir data/human/expression/latent_space_100/
fi

cd data/human/expression/latent_space_100/
human_latent_spaces=(../../../MLP_outcomes/*humantransform_1[0-9][0-9].csv)
for file in ${human_latent_spaces[@]}; 
do
    filename=$(basename -- "$file")
    if [ ! -f $filename ];
    then
        ln -s $file ./
    fi
done

#Return to main directory
cd ../../../..
