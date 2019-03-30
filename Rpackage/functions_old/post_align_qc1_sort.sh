#!/bin/bash

genome=$1
gtf=$2
bamin=$3
sample=$4
bamdir=$(dirname $bamin)

mkdir -p ${bamdir}/CollectMultipleMetrics
mkdir -p ${bamdir}/featureCounts

set -e

frag_dist=${bamdir}/CollectMultipleMetrics/${sample}.insert_size_metrics
gene_counts=${bamdir}/featureCounts/${sample} 

# check if fragent distribution and gene count files already exixt 
if [ -s $frag_dist -a -s $gene_counts ]; then 

	echo "$(basename $frag_dist) and $(basename $gene_counts) exist and are bigger than zero"

else

	sambamba sort ${bamin} -t $(nproc)

	bamin=${bamin/bam/sorted.bam}
	
	# collect multiple metrics
	echo CollectMultipleMetrics INPUT=${bamin} OUTPUT=${bamdir}/CollectMultipleMetrics/${sample} R=$genome 
	CollectMultipleMetrics INPUT=${bamin} OUTPUT=${bamdir}/CollectMultipleMetrics/${sample} R=$genome 

	# #gene counts

	echo featureCounts -p -T $(nproc) -a ${gtf} -o ${bamdir}/featureCounts/${sample} ${bamin}
	featureCounts -p -T $(nproc) -a ${gtf} -o ${bamdir}/featureCounts/${sample} ${bamin}

fi
