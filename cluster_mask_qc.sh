#!/bin/bash

species=mouse
#species=human

jacobians=absolute
#jacobians=relative

threshold=intensity
#threshold=topn

nk=3
k=1

if [ $jacobians == 'absolute' ];
then
	jacobians_short=abs
else
	jacobians_short=rel
fi

if [ $species == 'mouse' ]
then
  res=200
	template=data/mouse/atlas/DSURQE_CCFv3_average_200um.mnc
	map=data/${species}/clustering/cluster_maps/${jacobians}/resolution_${res}/mean/Group_${k}_Clusternum_${nk}_ES_${jacobians_short}_${res}_mean.mnc
else 
	template=data/human/registration/reference_files/model_1.0mm.mnc
	res=1.0
  map=data/${species}/clustering/cluster_maps/${jacobians}/resolution_${res}/mean/Group_${k}_Clusternum_${nk}_ES_${jacobians}_${res}_mean.mnc
fi

if [ $threshold == 'intensity' ];
then
	register $template \
	  $map \
		data/${species}/clustering/cluster_masks/${jacobians}/resolution_200/mean/threshold_${threshold}/symmetric/threshold_0.3/Group_${k}_Clusternum_${nk}_ES_${jacobians_short}_200_mean_mask_${threshold}_symmetric_signed_threshold0.3.mnc \
		data/${species}/clustering/cluster_masks/${jacobians}/resolution_200/mean/threshold_${threshold}/symmetric/threshold_0.5/Group_${k}_Clusternum_${nk}_ES_${jacobians_short}_200_mean_mask_${threshold}_symmetric_signed_threshold0.5.mnc \
		data/${species}/clustering/cluster_masks/${jacobians}/resolution_200/mean/threshold_${threshold}/symmetric/threshold_0.7/Group_${k}_Clusternum_${nk}_ES_${jacobians_short}_200_mean_mask_${threshold}_symmetric_signed_threshold0.7.mnc \
		data/${species}/clustering/cluster_masks/${jacobians}/resolution_200/mean/threshold_${threshold}/symmetric/threshold_0.9/Group_${k}_Clusternum_${nk}_ES_${jacobians_short}_200_mean_mask_${threshold}_symmetric_signed_threshold0.9.mnc
else
	register $template \
	  $map \
		data/${species}/clustering/cluster_masks/${jacobians}/resolution_200/mean/threshold_${threshold}/symmetric/threshold_0.2/Group_${k}_Clusternum_${nk}_ES_${jacobians_short}_200_mean_mask_${threshold}_symmetric_signed_threshold0.2.mnc \
		data/${species}/clustering/cluster_masks/${jacobians}/resolution_200/mean/threshold_${threshold}/symmetric/threshold_0.1/Group_${k}_Clusternum_${nk}_ES_${jacobians_short}_200_mean_mask_${threshold}_symmetric_signed_threshold0.1.mnc \
		data/${species}/clustering/cluster_masks/${jacobians}/resolution_200/mean/threshold_${threshold}/symmetric/threshold_0.05/Group_${k}_Clusternum_${nk}_ES_${jacobians_short}_200_mean_mask_${threshold}_symmetric_signed_threshold0.05.mnc

fi


