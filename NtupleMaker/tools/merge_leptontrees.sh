#!/bin/bash

# user input
tag=tnp_V00-00-00
path=/nfs-7/userdata/rwkelley/crab/lepton_trees 
merged_path=/nfs-7/userdata/rwkelley/babies/lepton_trees/$tag

# function to merge
function do_merge
{
	local input_path=$1
	local output_file=$2
	local treename="leptons"
	cmd="merge_tchain --input \"$input_path\" --output $output_file --tree $treename"
	echo $cmd
   	eval $cmd
}

# path to unmerged datasets 
unmerged_dataset_path=$path/$tag 

# get list of datasets
for dataset_dir in `ls -d $unmerged_dataset_path/*/`; 
do
	dataset=`basename ${dataset_dir}`
	log_name=merge_${dataset}.log
	mkdir -p logs
	do_merge "${dataset_dir}/res/leptonTree*.root" "${merged_path}/${dataset}/merged_leptontree.root" >& logs/${log_name} &
done
