#!/bin/bash

#Directory containing mouse cluster maps
datadir=/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters/

# Create symlinks for absolute 200um mean cluster maps
outdir=data/mouse/clustering/cluster_maps/absolute/resolution_200/mean
if [ ! -d '$outdir' ];
then 
    mkdir -p $outdir
fi

files_abs_200_mean=(${datadir}Group_*_ES_abs_200_mean.mnc)
for file in ${files_abs_200_mean[@]};
do
    ln -s $file ${outdir}
done


# Create symlinks for absolute 200um median cluster maps
outdir=data/mouse/clustering/cluster_maps/absolute/resolution_200/median
if [ ! -d '$outdir' ];
then 
    mkdir -p $outdir
fi

files_abs_200_median=(${datadir}Group_*_ES_abs_200_median.mnc)
for file in ${files_abs_200_median[@]};
do
    ln -s $file ${outdir}
done

# Create symlinks for relative 200um mean cluster maps
outdir=data/mouse/clustering/cluster_maps/relative/resolution_200/mean
if [ ! -d '$outdir' ];
then 
    mkdir -p $outdir
fi

files_rel_200_mean=(${datadir}Group_*_ES_rel_200_mean.mnc)
for file in ${files_rel_200_mean[@]};
do
    ln -s $file ${outdir}
done

# Create symlinks for relative 200um median cluster maps
outdir=data/mouse/clustering/cluster_maps/relative/resolution_200/median
if [ ! -d '$outdir' ];
then 
    mkdir -p $outdir
fi

files_rel_200_median=(${datadir}Group_*_ES_rel_200_median.mnc)
for file in ${files_rel_200_median[@]};
do
    ln -s $file ${outdir}
done

# Create symlink for cluster assignment .csv file
cluster_file=(${datadir}Clusters.csv)
ln -s $cluster_file data/mouse/clustering/mouse_clusters_groups10_200um.csv
