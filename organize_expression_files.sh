#!/bin/bash

if [ ! -d 'data/mouse/expression/latent_space_100/' ];
then
    mkdir data/mouse/expression/latent_space_100
fi

#Create symlinks for mouse latent space data
cd data/mouse/expression/latent_space_100/
mouse_latent_spaces=(../latent_space/*mousetransform_1[0-9][0-9].csv)
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

#Create human latent space directory
if [ ! -d "data/human/expression/latent_space_100/" ]; 
then
    mkdir data/human/expression/latent_space_100/
fi

cd data/human/expression/latent_space_100/
human_latent_spaces=(../latent_space/*humantransform_1[0-9][0-9].csv)
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
