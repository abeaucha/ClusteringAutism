#!/bin/bash

#Create mouse clustering directory
if [ ! -d 'data/mouse/clustering/' ]; 
then 
    mkdir data/mouse/clustering
fi 

#Create mouse cluster maps directory
if [ ! -d 'data/mouse/clustering/cluster_maps/' ];
then
    mkdir data/mouse/clustering/cluster_maps
fi

#Create mouse cluster maps absolute
if [ ! -d 'data/mouse/clustering/cluster_maps/absolute/' ];
then
    mkdir data/mouse/clustering/cluster_maps/absolute
fi

#Create mouse cluster maps absolute 200um dir
if [ ! -d 'data/mouse/clustering/cluster_maps/absolute/resolution_200' ];
then
    mkdir data/mouse/clustering/cluster_maps/absolute/resolution_200
fi

#Create mouse cluster maps absolute 200um mean dir
if [ ! -d 'data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/' ];
then
    mkdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean
fi

#Create symlinks for absolute 200um mean images
outdir=data/mouse/clustering/cluster_maps/absolute/resolution_200/mean
files_abs_200_mean=(/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters/Group_*_ES_abs_200_mean.mnc)

for file in ${files_abs_200_mean[@]};
do
    ln -s $file ${outdir}
done

#Create mouse cluster maps absolute 200um median dir
if [ ! -d 'data/mouse/clustering/cluster_maps/absolute/resolution_200/median/' ];
then
    mkdir data/mouse/clustering/cluster_maps/absolute/resolution_200/median
fi

#Create symlinks for absolute 200um median images
outdir=data/mouse/clustering/cluster_maps/absolute/resolution_200/median
files_abs_200_median=(/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters/Group_*_ES_abs_200_median.mnc)

for file in ${files_abs_200_median[@]};
do
    ln -s $file ${outdir}
done

#Create mouse cluster maps relative
if [ ! -d 'data/mouse/clustering/cluster_maps/relative/' ];
then
    mkdir data/mouse/clustering/cluster_maps/relative
fi

#Create mouse cluster maps relative 200um dir
if [ ! -d 'data/mouse/clustering/cluster_maps/relative/resolution_200' ];
then
    mkdir data/mouse/clustering/cluster_maps/relative/resolution_200
fi

#Create mouse cluster maps relative 200um mean dir
if [ ! -d 'data/mouse/clustering/cluster_maps/relative/resolution_200/mean/' ];
then
    mkdir data/mouse/clustering/cluster_maps/relative/resolution_200/mean
fi

#Create symlinks for relative 200um mean images
outdir=data/mouse/clustering/cluster_maps/relative/resolution_200/mean
files_rel_200_mean=(/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters/Group_*_ES_rel_200_mean.mnc)

for file in ${files_rel_200_mean[@]};
do
    ln -s $file ${outdir}
done

#Create mouse cluster maps relative 200um median dir
if [ ! -d 'data/mouse/clustering/cluster_maps/relative/resolution_200/median/' ];
then
    mkdir data/mouse/clustering/cluster_maps/relative/resolution_200/median
fi

#Create symlinks for relative 200um median images
outdir=data/mouse/clustering/cluster_maps/relative/resolution_200/median
files_rel_200_median=(/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters/Group_*_ES_rel_200_median.mnc)

for file in ${files_rel_200_median[@]};
do
    ln -s $file ${outdir}
done

cluster_file=/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters/Clusters.csv

ln -s $cluster_file data/mouse/clustering/mouse_clusters_groups10_200um.csv

