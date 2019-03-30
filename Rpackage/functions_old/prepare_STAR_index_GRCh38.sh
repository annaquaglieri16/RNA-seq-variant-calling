#!/bin/bash

##########################################
## Prepare files for star fusion with hg38
##########################################
workdir=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/quaglieri.a/scripts/functions
genomedir=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/
genome_fasta_ucsc=${genomedir}/hg38/UCSC/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
genome_fasta_ncbi=${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
gtf=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/UCSC/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf
gtf_ncbi=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
star_fusion_data_ucsc=${genomedir}/star_fusion_datasets_UCSC
star_fusion_data_ncbi=${genomedir}/star_fusion_datasets_NCBI
# remove the unknown chromosomes since they are not present in the genome.fa
# egrep "^chr[0-9XYM]{1,2}\s" ${gtf} > /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/UCSC/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes_known_chrom.gtf
# egrep "^chr[0-9XYM]{1,2}\s" ${gtf_ncbi} > /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes_known_chrom.gtf
gtf_known=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/UCSC/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes_known_chrom.gtf
gtf_known_ncbi=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes_known_chrom.gtf

# UCSC reference fasta
${workdir}/star_fusion_prepare.sh ${gtf_known} ${genome_fasta_ucsc} ${star_fusion_data_ucsc}

# NCBI reference fasta - unix401
${workdir}/star_fusion_prepare.sh ${gtf_known_ncbi} ${genome_fasta_ncbi} ${star_fusion_data_ncbi}


###################################################
## Build STAR index for the 100bp libraries GRChg38
###################################################

readlen_1=99
genomedir=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes
genome_fasta_ucsc=${genomedir}/hg38/UCSC/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
genome_fasta_ncbi=${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
gtf=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/UCSC/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf
gtf_ncbi=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
workdir=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/quaglieri.a/scripts/functions

# UCSC fasta
assembly="hg38"
${workdir}/build_STAR_index_P1.sh ${genome_fasta_ucsc} ${gtf} ${genomedir} ${assembly} ${readlen_1}

# NCBI fasta
assembly="GRCh38"
${workdir}/build_STAR_index_P1.sh ${genome_fasta_ncbi} ${gtf_ncbi} ${genomedir} ${assembly} ${readlen_1}

