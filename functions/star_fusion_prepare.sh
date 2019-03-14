#!/bin/bash

# ## Building a custom FUsionFilter Dataset
# # https://github.com/FusionFilter/FusionFilter/wiki/Building-a-Custom-FusionFilter-Dataset
# Depends on the Trinity software
# Downalod STAR-Fusion from GitHub

# datasets
#gtf=$1
#genome_fasta=$2
#star_fusion_data=$3
#star_fusion_source=$4

star_fusion_source=/wehisan/home/allstaff/q/quaglieri.a/software/STAR-Fusion-v1.4.0
gtf=/wehisan/home/allstaff/q/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes_known_chrom.gtf
ve-2015-08-11-09-31-31/Genes/genes.gtf
genome_fasta=/wehisan/home/allstaff/q/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

mkdir -p ${star_fusion_data}

cd ${star_fusion_data}

# 0. Download fusion annotation
# check if fusion annotation is present or get it from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
# annotation=./Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v26_CTAT_lib_July192017.source_data/fusion_lib.dat.gz
annotation=/wehisan/home/allstaff/q/quaglieri.a/software/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v27_CTAT_lib_July192017.source_data/fusion_lib.dat.gz

if [[ -r ${annotation} ]]; then 

	echo "Fusion annotation has already been downloaded"

else


  wget -r -np --no-host-directories https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.source_data.tar.gz -P .

	#wget -r -np --no-host-directories https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v26_CTAT_lib_July192017.source_data.tar.gz -P .

	tar xvzf /wehisan/home/allstaff/q/quaglieri.a/software/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018/GRCh38_v27_CTAT_lib_Feb092018.source_data.tar.gz 

	annotation=/wehisan/home/allstaff/q/quaglieri.a/software/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018/CTAT_HumanFusionLib.v0.1.0.dat.gz

fi


## Extract all _alt features in gtf files that are not in the fasta file


# New way of preparing the files fro STAR-Fusion
ctat_dir=/wehisan/home/allstaff/q/quaglieri.a/software/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018/
cd $ctat_dir

module load STAR/2.5.3a
module load trinity/2.4.0
module load samtools/1.9

$star_fusion_source/FusionFilter/prep_genome_lib.pl \
                         --genome_fa $genome_fasta \
                         --gtf $gtf \
                         --fusion_annot_lib CTAT_HumanFusionLib.dat.gz \
                         --annot_filter_rule AnnotFilterRule.pm \
                         --pfam_db PFAM.domtblout.dat.gz
