#!/bin/bash

seqtkdir=/home/users/allstaff/quaglieri.a/software/seqtk
fastqin=$1
sample=$2
fastqdown_dir=$3
nreads=$4
seed=$5

set -e

mkdir -p ${fastqdown_dir}

if [[ -s ${fastqdown_dir}/${sample}_R1.fastq.gz ]] && [[ -s ${fastqdown_dir}/${sample}_R2.fastq.gz ]] ; then

	echo "${fastqdown_dir}/${sample}_R1.fastq and ${fastqdown_dir}/${sample}_R2.fastq already exixts"

else

	echo "downsampling $(basename $fastqin)"

	FQ1=${fastqin}
	echo "$FQ1"
	FQ2=${FQ1/_1.fastq/_2.fastq}
	echo "$FQ2"

	${seqtkdir}/seqtk sample -s${seed} $FQ1 ${nreads} > ${fastqdown_dir}/${sample}_R1.fastq

	${seqtkdir}/seqtk sample -s${seed} $FQ2 ${nreads} > ${fastqdown_dir}/${sample}_R2.fastq

	#pigz -p $(nproc) -4 ${fastqdown_dir}/${sample}_R1.fastq
	#pigz -p $(nproc) -4 ${fastqdown_dir}/${sample}_R2.fastq

	gzip ${fastqdown_dir}/${sample}_R1.fastq
	gzip ${fastqdown_dir}/${sample}_R2.fastq


fi



