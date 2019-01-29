#!/bin/bash

module load STAR

star_fusion=/home/users/allstaff/quaglieri.a/software/STAR-Fusion
genomedir=$1
gtf=$2
chim_junction_dir=$3
sample=$4
outdir=$5

mkdir -p ${outdir}

out=${outdir}/${sample}

if [[ -s ${out}/star-fusion.fusion_predictions.abridged.tsv ]] ; then 

	echo "$out already exists"

else

	${star_fusion}/STAR-Fusion \
	--genome_lib_dir ${genomedir} \
	-J ${chim_junction_dir}/${sample}Chimeric.out.junction \
	--CPU $(nproc) \
	--output_dir ${outdir}/${sample}

fi