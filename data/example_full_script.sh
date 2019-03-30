#!/bin/bash

################
## Load Modules
################

module load sambamba
module load gmap-gsnap/2016-04-04
module load samtools/1.3.1
module load gatk/3.7.0


#########################################################################
##### Directories and file names: NCBI fasta file and NCBI annotation GTF
#########################################################################

# downloaded all UCSC hg38 from iGenomes 
# website https://support.illumina.com/sequencing/sequencing_software/igenome.html

star_index="star_index_GRCh38_99"
star_fusion="star_fusion_datasets_NCBI"
assembly="GRCh38"
# Genome directories
genomedir=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes
genome_fasta=${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
gtf=${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
star_genome100=${genomedir}/${star_index}
star_fusion_data_ncbi=${genomedir}/${star_fusion}

### Working directories
workdir=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/quaglieri.a/scripts/functions
fastqdir=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/fastq
amldir=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf
star_1pass=${amldir}/aligned_pass1_${assembly}
star_2pass=${amldir}/aligned_pass2_${assembly}
star_merged=${star_2pass}_merged_runs

# Fastq files
samplein_runs=$(find ${fastqdir} -maxdepth 1 -name "*R1.fastq*" | sort)

# create files with names for output files
runs=$(find ${fastqdir} -maxdepth 1  -name "*fastq.gz" | cut -f9 -d '/' | cut -f1,2,4 -d "_" | sort | uniq)
samples=$(find ${fastqdir} -maxdepth 1  -name "*fastq.gz" | cut -f9 -d '/' | cut -f1 -d "_" | sort | uniq)
find ${fastqdir} -maxdepth 1  -name "*fastq.gz" | cut -f9 -d '/' | cut -f1,2,4 -d "_" | sort | uniq > ${amldir}/runs.txt
cat ${amldir}/runs.txt | cut -f1 -d "_" | sort | uniq  > ${amldir}/samples_unique.txt
cat ${amldir}/runs.txt | cut -f1 -d "_"  > ${amldir}/samples.txt
paste ${amldir}/samples.txt ${amldir}/runs.txt > ${amldir}/samples_runs.txt
paste ${amldir}/samples.txt ${amldir}/samples.txt > ${amldir}/samples_samples.txt
runs=${amldir}/runs.txt
samples_runs=${amldir}/samples_runs.txt
# SRX729626 SRR1608879
# SRX729626 SRR1608881
# SRX729626 SRR1608885
# SRX729622 SRR1608851
# SRX729622 SRR1608850
samples_samples=${amldir}/samples_samples.txt

# GATK and variant calling required files 
mills1k_indels=${genomedir}/hg38/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
dbSNP=${genomedir}/hg38/gatk_bundle/dbsnp_144.hg38_chr.vcf.gz
snps1k=${genomedir}/hg38/gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz
regions=${amldir}/mutated_genes_aml/genes_for_variant_calling.bed
outvariants=${amldir}/variant_calling

## HPC script_dir
script_dir=${fastqdir}/../HPC_scripts


###################################
## Fastqc analysis - Quality checks
find ${fastqdir} -name "*.fastq*" > ${fastqdir}/fastq_files.txt
cat ${fastqdir}/fastq_files.txt | parallel -j 10 "fastqc {} --outdir ${fastqdir}/fastqc"

## count how many zip files have been created
# find ${fastqdir}/fastqc -maxdepth 1 -name "*.zip" | wc -l
# #  zipped files

# Multi QC report
multiqc ${fastqdir}/fastqc --interactive -n "mydata_summary_report" -o ${fastqdir}/../ 

##############
##############
# STAR-pass1 
##############
##############

for fastqin in ${samplein_runs} ; do

# Define sample name for outpur file
bamout=$(echo $(basename $fastqin) | cut -f1 -d "_")
# bamout=${fastqin/_R1.fastq.gz/}

echo '#!/bin/bash' > ${script_dir}/${bamout}_star1.sh
echo 'module load STAR' >> ${script_dir}/${bamout}_star1.sh
echo '' >> ${script_dir}/${bamout}_star1.sh

# STAR input fastq files
FQ1=$fastqin
FQ2=${FQ1/R1/R2}

echo ${workdir}/run_STAR_align_P1.sh \
${star_genome100} \
${FQ1} ${FQ2} \
${bamout} \
${star_1pass} >> ${script_dir}/${bamout}_star1.sh

echo chmod 775 ${star_1pass} >> ${script_dir}/${bamout}_star1.sh

done

###########
# On torque
fastqdir=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/fastq
script_dir=${fastqdir}/../HPC_scripts

scripts=$(find ${script_dir} -name "*star1.sh")

for script in $scripts; do

	cd ${script_dir}
	qsub -n -q submit -l nodes=1:ppn=20,mem=35gb ${script}

done

###########

## Check alignment
# Multi QC report
multiqc ${fastqdir}/fastqc ${star_1pass} --interactive -n "mydata_summary_report" -f -o ${star_1pass}/..


# concatenate splice junctions from all samples from ran in pass1
cat ${star_1pass}/*SJ.out.tab > ${star_1pass}/combined_sj.out.tab
# Dobin suggests to remove chrm cause they are usually False positives
awk '!/chrM/' ${star_1pass}/combined_sj.out.tab > ${star_1pass}/combined_sj_nochrM.out.tab

######################################################
# Check integrity of bamfiles by looking at the header of the aligned bamfile
bams=$(find ${star_1pass} -maxdepth 1 -name "*Aligned.sortedByCoord.out.bam" | sort)
for bam in $bams ; do
echo $bam
samtools view $bam -H | head -1
done

#############
## STAR pass2
#############

for fastqin in ${samplein_runs} ; do

# Define sample name for outpur file
bamout=$(echo $(basename $fastqin) | cut -f1 -d "_")
# bamout=${fastqin/_R1.fastq.gz/}
# sampleID=$(echo $(basename $fastqin) | cut -f1 -d "_") if more fastq files are going to be merged together

echo '#!/bin/bash' > ${script_dir}/${bamout}_star2.sh
echo 'module load STAR' >> ${script_dir}/${bamout}_star2.sh
echo '' >> ${script_dir}/${bamout}_star2.sh

# STAR input fastq files
FQ1=$fastqin
FQ2=${FQ1/R1/R2}

echo ${workdir}/run_STAR_align_P2.sh \
${star_genome100} \
$FQ1 $FQ2 \
${bamout} \
${star_2pass} \
${star_1pass}/combined_sj_nochrM.out.tab >> ${script_dir}/${bamout}_star2.sh
echo '' >> ${script_dir}/${bamout}_star2.sh

######################################
# Mark duplicates, addRG and Validate 
######################################
bamin=${star_2pass}/${bamout}Aligned.sortedByCoord.out.bam

echo ${workdir}/post_align_qc1.sh ${genome_fasta} ${gtf} ${star_2pass} ${bamout} >> ${script_dir}/${bamout}_star2.sh
echo '' >> ${script_dir}/${bamout}_star2.sh

echo ${workdir}/post_align_qc2.sh ${bamin} ${bamout} ${genome_fasta} ${bamout} >> ${script_dir}/${bamout}_star2.sh

done

######################################################
######################################################
# Check integrity of bamfiles by looking at the header of the aligned bamfile
find . -maxdepth 1 -name "*_Aligned.reorderedDupl.rg.split.bam" | wc -l 
bams=$(find ${star_2pass} -maxdepth 1 -name "*_Aligned.reorderedDupl.rg.split.bam" | sort)
for bam in $bams ; do
echo $bam
samtools view $bam -H | head -1
done

## Summary of validation
Rscript ${workdir}/readValidateSamFile.R ${star_2pass}/ValidateSam

# Combine fragment sizes and PCA of samples
Rscript ${workdir}/combine_frag_MDSplot.R \
${star_2pass}/CollectMultipleMetrics ${star_2pass}/featureCounts ${samples_runs} "all_runs"

# Multi QC report
multiqc ${fastqdir}/fastqc ${star_2pass} --interactive -f -n "mydata_summary_report" -o ${star_2pass}/../ 


###################################################################
# Merge runs & Combine chimeric reads from runs of the same sample 
###################################################################

for sample in ${samples} ; do

sample=${sample}_

echo '#!/bin/bash' > ${script_dir}/${sample}merge_and_fusion.sh
echo 'module add samtools' >> ${script_dir}/${sample}merge_and_fusion.sh
echo 'module add sambamba' >> ${script_dir}/${sample}merge_and_fusion.sh
echo 'module add STAR' >> ${script_dir}/${sample}merge_and_fusion.sh
echo 'module add R' >> ${script_dir}/${sample}merge_and_fusion.sh
echo '' >> ${script_dir}/${sample}merge_and_fusion.sh

echo ${workdir}/merge_runs.sh ${sample} ${star_2pass} >> ${script_dir}/${sample}merge_and_fusion.sh

# Fusion detection
# merge chimeric reads from runs of the same sample
cat ${star_2pass}/${sample}*L00*Chimeric.out.junction > ${script_dir}/${sample}Chimeric.out.junction

# Run STAR-Fusion
echo ${workdir}/star_fusion.sh \
${star_fusion_data_ncbi} ${gtf} ${star_2pass} ${sample} >> ${script_dir}/${sample}merge_and_fusion.sh

echo $sample

done

# Run on torquewelord
ssh torquelord1.hpc.wehi.edu.au

scripts=$(find ${script_dir} -maxdepth 1 -name '*_merge_and_fusion.sh')

for script in $scripts ; do

	cd ${script_dir}
	qsub -n -q submit -l nodes=1:ppn=20,mem=35gb ${script}

done

########################################################
# Check integrity of merge files through count of proper pairs in flagstat
########################################################

Rscript --vanilla ${workdir}/flagstat_summary.R --directory ${star_2pass} --pattern "Aligned\\.reorderedDupl\\.rg\\.bam$"

Rscript --vanilla ${workdir}/flagstat_summary.R --directory ${star_merged} --pattern "_merged_Dupl_rg\\.bam$"

R
comb_before <- read.csv("/wehisan/general/user_managed/grpu_majewski_3/venetoclax_trial/GRCh38_aligned/aligned_pass2/flagstats/combine_flagstats.csv")
comb_merged <- read.csv("/wehisan/general/user_managed/grpu_majewski_3/venetoclax_trial/GRCh38_aligned/aligned_pass2_merged_runs/flagstats/combine_flagstats.csv")

library(dplyr)
library(tidyr)
comb_before_sum <- comb_before %>% tidyr::separate(ID,into=c("ID","barcode","lane")) %>%
group_by(ID) %>% dplyr::summarise(libsizeBefore = sum(properPaired))

combine  <- merge(comb_merged,comb_before_sum)
summary(combine$diff)

q()

###############
# plot fusions
###############
R
data <- read.table("/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/samples_runs.txt",stringsAsFactors=FALSE)
data[,1] <- paste(data[,1],"_",sep="")
write.table(data,"/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/samples_runs1.txt",row.names=F,col.names=F,sep="\t",quote=FALSE)
q()
n

samples_runs1=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/samples_runs1.txt

Rscript ${workdir}/analyse_fusion_output2.R \
${star_2pass}/star_fusion $samples_runs1 ${star_2pass}/featureCounts 5

Rscript ${workdir}/analyse_fusion_output2.R \
${star_2pass}/star_fusion $samples_runs1 ${star_2pass}/featureCounts 10


#################################
#####  Bams Pre-Process ##########
#################################
## GATK pre-processes 

# Input needs to have RG and Duplicates marked
bams=$(find ${star_merged} -maxdepth 1 -name "*_merged_Dupl_rg.bam" | sort)

for bamin in $bams; do

sample_tmp=$(basename $bamin)
sample=${sample_tmp/_merged_Dupl_rg.bam/}

echo '#!/bin/bash' > ${script_dir}/${sample}_gatk_pipe.sh
echo 'module add gatk' >> ${script_dir}/${sample}_gatk_pipe.sh
echo 'module add R' >> ${script_dir}/${sample}_gatk_pipe.sh
echo 'module add sambamba' >> ${script_dir}/${sample}_gatk_pipe.sh
echo 'module add varscan' >> ${script_dir}/${sample}_gatk_pipe.sh
echo 'module add vcftools' >> ${script_dir}/${sample}_gatk_pipe.sh
echo 'module add samtools' >> ${script_dir}/${sample}_gatk_pipe.sh
echo 'module add ensembl-vep' >> ${script_dir}/${sample}_gatk_pipe.sh
echo 'set -e' >> ${script_dir}/${sample}_gatk_pipe.sh
echo ' ' >> ${script_dir}/${sample}_gatk_pipe.sh
echo '# GATK pre-process BAM ' >> ${script_dir}/${sample}_gatk_pipe.sh


######################
## GATK pre-processes 
######################

echo Rscript ${workdir}/gatk_process_pipe.R --reference_fasta ${genome_fasta} \
--bamfile ${bamin}  --sampleName ${sample} --knownSites2 ${dbSNP} --knownSites1 ${indels1k} \
--knownSites3 ${mills1k_indels} --run TRUE >> ${script_dir}/${sample}_gatk_pipe.sh

bamin_final=${bamin/Aligned.reorderedDupl.rg.bam/Recal.reorderedDupl.rg.split.bam}

############################
# call and annotate variants
############################

echo ' ' >> ${script_dir}/${sample}_gatk_pipe.sh
echo '# MUtect call variant ' >> ${script_dir}/${sample}_gatk_pipe.sh
echo ' ' >> ${script_dir}/${sample}_gatk_pipe.sh

echo Rscript ${workdir}/call_variants.R --reference_fasta ${genome_fasta} --genome_assembly "GRCh38" \
--bamfile ${bamin_final} --sampleName ${sample} --dbSNP ${dbSNP} -c "mutect" --output_directory ${outvariants} \
--VEPcache "/stornext/HPCScratch/cache/.vep" >> ${script_dir}/${sample}_gatk_pipe.sh

echo ' ' >> ${script_dir}/${sample}_gatk_pipe.sh
echo '# VarScan call variant ' >> ${script_dir}/${sample}_gatk_pipe.sh
echo ' ' >> ${script_dir}/${sample}_gatk_pipe.sh

echo Rscript ${workdir}/call_variants.R --reference_fasta ${genome_fasta} --genome_assembly "GRCh38" \
--bamfile ${bamin_final} --sampleName ${sample} -c "varscan" --output_directory ${outvariants} \
--VEPcache "/stornext/HPCScratch/cache/.vep" >> ${script_dir}/${sample}_gatk_pipe.sh

done


# Run on toriwelord
ssh torquelord4.hpc.wehi.edu.au

scripts=$(find ${script_dir} -maxdepth 1 -name '*gatk_pipe.sh')
cd ${script_dir}

for script in $scripts ; do

cd ${script_dir}
qsub -n -q submit -l nodes=1:ppn=3,mem=64gb ${script}

done


