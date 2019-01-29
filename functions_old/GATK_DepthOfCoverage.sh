#!/bin/bash

# Create scripts to compute coverage for every downsampled run 

run=$1
bam_find=$2
loci=$3
# /wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/scripts/02-variant_loss_downsampling/mutations_loci.bed

bamdir=/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/hg19_aligned/STAR_2pass_${run}_AML
find $bamdir -name "*${bam_find}*bam" > ${bamdir}/list_bams_${run}.list
script_dir=/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/scripts/06-Streamline_downsampling/05-DepthOfCoverage
mkdir -p ${script_dir}/${run}

loci=$3
genome_fasta=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
bams=${bamdir}/list_bams_${run}.list

echo '#!/bin/bash' > ${script_dir}/${run}_doc.sh
echo ' ' >> ${script_dir}/${run}_doc.sh

echo '#PBS -q submit' >> ${script_dir}/${run}_doc.sh
echo '#PBS -l nodes=1:ppn=2,mem=16gb' >> ${script_dir}/${run}_doc.sh
echo '#PBS -l walltime=00:10:00' >> ${script_dir}/${run}_doc.sh
echo '#PBS -o ' ${script_dir}/${run}_doc_out >> ${script_dir}/${run}_doc.sh
echo '#PBS -e ' ${script_dir}/${run}_doc_err >> ${script_dir}/${run}_doc.sh
#echo '#PBS -N germCall_CBF-AML_WEHI' >> ${script_dir}/${run}_doc.sh
echo ' ' >> ${script_dir}/${run}_doc.sh


echo 'module add gatk/3.7.0' >> ${script_dir}/${run}_doc.sh
echo 'set -e' >> ${script_dir}/${run}_doc.sh
echo ' ' >> ${script_dir}/${run}_doc.sh

echo gatk -T DepthOfCoverage \
-R ${genome_fasta} \
-o ${script_dir}/${run}/${run} \
-I ${bams} \
--printBaseCounts \
-L ${loci} >> ${script_dir}/${run}_doc.sh

