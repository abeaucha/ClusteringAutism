#!/bin/bash
# coordinates_to_tags.sh
# Author: Antoine Beauchamp
# Convert a CSV of coordinates to Tag point format

outdir=data/human/expression/
# coordinates=${outdir}AHBA_microarray_coordinates_mni.csv
# coordinates=${outdir}AHBA_microarray_coordinates_study_v1.csv
coordinates=${outdir}AHBA_microarray_coordinates_study_v2.csv

tags=$(basename $coordinates .csv).tag
tags=${outdir}$tags

echo """MNI Tag Point File
Volumes = 1;

Points = """ > $tags

tail -n +2 $coordinates | awk -v FPAT="([^,]+)|(\"[^\"]+\")" -v OFS=" " '{print $1,$2,$3,1,$5,1,$5}' >> $tags

echo ";" >> $tags
