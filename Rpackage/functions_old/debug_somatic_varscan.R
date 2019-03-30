# Run variant calling 

# Get Counts with FeatureCounts 
local <- FALSE
dir <- ifelse(local,"/Volumes/AML_PROJECT/AML_RNA/cbf_aml_agrf","/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf")
mydir <- ifelse(local,"/Volumes/Anna's\ UNIX\ Home\ area/PHD_project/GEO_Leucegene_data","/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data")

# GATK and variant calling required files 
workdir <- file.path(dir, "../quaglieri.a/scripts/functions")
star_merged <- file.path(dir,"aligned_pass2_GRCh38_merged_runs")
genome_fasta <- file.path(mydir,"genomes","hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa")

mills1k_indels <- file.path(mydir,"genomes","hg38/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
dbSNP <- file.path(mydir,"genomes","hg38/gatk_bundle/dbsnp_144.hg38.vcf.gz")
snps1k <- file.path(mydir,"genomes","hg38/gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz")
regions <- file.path(dir,"mutated_genes_aml/genes_for_variant_calling.bed")
outvariants <- file.path(dir,"variant_calling")
normal_vcf <- file.path(mydir,"hg19_aligned/variant_calling/samtools/combined_variants_control_cd341_filtered.vcf")


# Get Matched table
matched_table <- read.table(file.path(dir,"variant_calling","cancer_normal_matches"),stringsAsFactors = FALSE,header = TRUE)
i=1
# Get list of bamfiles
# bamfiles <- list.files(path=star_merged,pattern="_Recal\\.reorderedDupl\\.rg\\.split\\.bam$")
# dim(matched_table) # 42 including replicates as samples in diagnosis and relapse from Data Structure

for( i in 1:nrow(matched_table)) {

	cancer <- as.character(matched_table[i,1])
	matched <- as.character(matched_table[i,2])
	normal <- ifelse(is.na(matched),NULL,matched)

	sample <- gsub("_Recal.reorderedDupl.rg.split.bam","",cancer)

	# Create lock for bamfile
	lock <- file.path(star_merged,paste0(sample,"_processed.lock"))

	if ( dir.exists(lock) ){

		print(paste0(cancer, " is being processed."))

	} else {

		dir.create(lock,recursive = TRUE,showWarnings = FALSE)

		system(paste0("Rscript ", file.path(workdir,"call_variants.R"),
			" --reference_fasta ", genome_fasta,
			" --bamfile ", file.path(star_merged,cancer),
			" --matched_normal ", file.path(star_merged,normal),
			" --sampleName ", sample,
			" --dbSNP ", dbSNP,
			" -c mutect",
			" --output_directory ", outvariants,
			" --regions ", regions))

		unlink(lock,recursive = FALSE)

	}

}


system(paste0("Rscript ", file.path(workdir,"call_variants.R"),
			" --reference_fasta ", genome_fasta,
			" --bamfile ", file.path(star_merged,cancer),
			" --matched_normal ", file.path(star_merged,normal),
			" --sampleName ", sample,
			" --dbSNP ", dbSNP,
			" -c varscan",
			" --output_directory ", outvariants,
			" --regions ", regions))


samtools mpileup --count-orphans --min-BQ 20 --min-MQ 20 -m 3 -F 0.0002 \
--fasta-ref /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
-l /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/mutated_genes_aml/genes_for_variant_calling.bed \
 /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38_merged_runs/4_Recal.reorderedDupl.rg.split.bam \
 /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38_merged_runs/2_Recal.reorderedDupl.rg.split.bam

--count-orphans --min-BQ 20 --min-MQ 20 --output-BP --output-MQ -m 3 -F 0.0002

samtools mpileup --count-orphans --min-BQ 20 --min-MQ 20 -m 3 -F 0.0002 --output-tags AD,ADF,ADR,DP,SP \
--fasta-ref /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
-l /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/mutated_genes_aml/genes_for_variant_calling.bed \
 /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38_merged_runs/4_Recal.reorderedDupl.rg.split.bam \
 /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38_merged_runs/2_Recal.reorderedDupl.rg.split.bam 


 | VarScan somatic \
-mpileup "somatic" --min-coverage 10 --min-reads2 2 --variants 1 --output-vcf 1 --min-avg-qual 20 --min-var-freq 0.01 


# vep

vcf=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/mutect/regions/9_somatic_snvs_indels.vcf


# without --output-BP
chr1    36465544        T       34      ,,,,,,,,.......,,..,..,,...,...,,.      IA@@@@@@@_^_`a_@@@@@@@@@__@@_@@@>=      33      ,,,,.,,,.,.,..,.,,,.,,...,.,.....       DBB
AABBBAAAB@AB`B@BABBAAABABA=>8=

# with --output-BP
chr1    36465544        T       34      ,,,,,,,,.......,,..,..,,...,...,,.      IA@@@@@@@_^_`a_@@@@@@@@@__@@_@@@>=      99,93,89,87,82,82,82,81,80,80,75,74,74,74,74,72,69,
69,68,68,68,64,64,64,64,62,62,59,58,55,55,45,40,38,38,37,35,35,35,29,20,13,3,2  33      ,,,,.,,,.,.,..,.,,,.,,...,.,.....       DBBAABBBAAAB@AB`B@BABBAAABABA=>8=       97,
91,82,80,79,77,71,71,68,64,64,62,62,59,56,56,51,41,37,37,37,32,32,31,30,27,21,21,19,19,13,2,2,3,2


###############
#### Debug VEP
###############

vcf=/home/users/allstaff/quaglieri.a/9_somatic_snvs_indels.vcf

# Annotate variants
# if (annotate) {

# 	vep_call <- paste0("vep --dir_cache /usr/local/bioinfsoftware/VEP/VEP-85/variant_effect_predictor.pl",
# 				"-i ", VCF, " -o ",annotatedVCF, " --cache --everything --force_overwrite ", "--fork 12")

# 	print(vep_call)
# 	system(vep_call)
# }

cp /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/mutated_genes_aml/genes_for_variant_calling.bed /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/mutated_genes_aml/genes_for_variant_calling_copy.bed
regions_copy=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/mutated_genes_aml/genes_for_variant_calling_copy.bed
mutated_genes=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/mutated_genes_aml
gtf=${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf

cp ${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf ${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes_copy.gtf

cd /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/mutated_genes_aml
bgzip -f genes_for_variant_calling_copy.bed
tabix -p bed genes_for_variant_calling_copy.bed.gz 

bgzip -f ${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes_copy.gtf
tabix -p gff ${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes_copy.gtf.gz
gtf_gz=${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes_copy.gtf.gz




##################
#####
genomedir=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes
genome_fasta=${genomedir}/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
var=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/varscan/regions/9_somatic_snvs_indels.vcf.snp

cp /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/varscan/regions/9_somatic_snvs_indels.vcf.snp \
/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/varscan/regions/9_somatic_snvs_indels.snp.vcf
var1=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/varscan/regions/9_somatic_snvs_indels.snp.vcf

head -40 $var > ~/subset_var.vcf
subset=~/subset_var.vcf
# VEP
/usr/local/bioinfsoftware/VEP/VEP-85/bin/variant_effect_predictor.pl --dir_cache /usr/local/bioinfsoftware/VEP/VEP-85/.vep \
--assembly GRCh38 -i ${subset} -o ~/annotatedVCF.vcf --cache --everything --force_overwrite --vcf 

/usr/local/bioinfsoftware/VEP/VEP-85/bin/variant_effect_predictor.pl --dir_cache /usr/local/bioinfsoftware/VEP/VEP-85/.vep \
--assembly GRCh38 -i ${subset} -o ~/annotatedVCF_tab --cache --everything --force_overwrite --tab

##### Varscan 
vcf-query ${var1} -f '%CHROM %POS %REF %ALT %FILTER %INFO/SPV %INFO/GPV [ %GT %GQ %DP %FREQ %RD %AD %DP4] \n' | head

# Mutect2
mut=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/mutect/regions/16_somatic_snvs_indels.vcf
vcf-query ${mut} -f '%CHROM %POS %REF %ALT %FILTER [ %GT %QSS %AD %AF %ALT_F1R2 %ALT_F2R1 %REF_F1R2 %REF_F2R1] \n' > ~/mutect_fields


mut_old=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/variant_calling/hg19/mutect2/SRX958862_snvs_indels.vcf
vcf-query ${mut_old} -f '%CHROM %POS %REF %ALT %FILTER [ %GT %QSS %AD %AF %ALT_F1R2 %ALT_F2R1 %REF_F1R2 %REF_F2R1] \n' > ~/mutect_old_fields

var_old=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/variant_calling/hg19/down_20M/samtools/SRX958862Recal_real_varscan_snvs.vcf
vcf-query ${var_old} -f '%CHROM %POS %REF %ALT %FILTER [ %GT %RBQ %ABQ %DP %FREQ %RD %AD %RDF %RDR %ADF %ADR] \n' > ~/var_old_fields

vcf_ours <- read.table("~/mutect_fields")[,1:13]
vcf <- read.table("~/mutect_old_fields")



##FILTER=<ID=alt_allele_in_normal,Description="Evidence seen in the normal sample">
##FILTER=<ID=clustered_events,Description="Clustered events observed in the tumor">
##FILTER=<ID=clustered_read_position,Description="Evidence for somatic variant clusters near the ends of reads">
##FILTER=<ID=germline_risk,Description="Evidence indicates this site is germline, not somatic">
##FILTER=<ID=homologous_mapping_event,Description="More than three events were observed in the tumor">
##FILTER=<ID=multi_event_alt_allele_in_normal,Description="Multiple events observed in tumor and normal">
##FILTER=<ID=panel_of_normals,Description="Seen in at least 2 samples in the panel of normals">
##FILTER=<ID=str_contraction,Description="Site filtered due to contraction of short tandem repeat region">
##FILTER=<ID=strand_artifact,Description="Evidence for alt allele comes from one read direction only">
##FILTER=<ID=t_lod_fstar,Description="Tumor does not meet likelihood threshold">
##FILTER=<ID=triallelic_site,Description="Site filtered because more than two alt alleles pass tumor LOD">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of the event in the tumor">
##FORMAT=<ID=ALT_F1R2,Number=1,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting the alternate allele">
##FORMAT=<ID=ALT_F2R1,Number=1,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting the alternate allele">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=FOXOG,Number=1,Type=Float,Description="Fraction of alt reads indicating OxoG error">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=QSS,Number=A,Type=Integer,Description="Sum of base quality scores for each allele">
##FORMAT=<ID=REF_F1R2,Number=1,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting the reference allele">
##FORMAT=<ID=REF_F2R1,Number=1,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting the reference allele">

######################################################
## Merge VCF so to have one VCF from the varscan output
indel=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/varscan/regions/9_somatic_snvs_indels.vcf.indel
cp $indel ${indel}.vcf
snp=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/varscan/regions/9_somatic_snvs_indels.vcf.snp
cp $snp ${snp}.vcf
dir=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/varscan/regions

vcf-concat ${indel}.vcf ${snp}.vcf > ${dir}/concat.vcf

/usr/local/bioinfsoftware/VEP/VEP-85/bin/variant_effect_predictor.pl --dir_cache /usr/local/bioinfsoftware/VEP/VEP-85/.vep \
--assembly GRCh38 -i ${dir}/concat.vcf -o ~/annotatedVCF_tab_concat --cache --everything --force_overwrite --tab


########################################
## Call variants on the whole genome
########################################

opt <- data.frame(reference_fasta=genome_fasta,bamfile=file.path(star_merged,cancer),
	matched_normal=as.character(file.path(star_merged,normal)),sampleName=sample,dbSNP=dbSNP,
	variant_caller="mutect",output_directory=as.character(outvariants),regions=regions)

opt <- data.frame(reference_fasta=genome_fasta,bamfile=file.path(star_merged,cancer),
	matched_normal=as.character(file.path(star_merged,normal)),sampleName=sample,dbSNP=dbSNP,
	variant_caller="varscan",output_directory=as.character(outvariants),regions=regions)

# call variants on the whole genome
opt <- data.frame(reference_fasta=genome_fasta,bamfile=file.path(star_2pass,bamfile),
	sampleName=sample,dbSNP=dbSNP,matched_normal=NA,
	variant_caller="mutect",output_directory=as.character(outvariants), create_PON=1)

mutect_call <- paste0("gatk -T MuTect2",
			" -R ", genome, 
			" --read_filter MappingQuality --min_mapping_quality_score 20 --max_alternate_alleles 4 -kmerSize 15 -kmerSize 30", 
			" -I:tumor ", bamfile)

mutect_call <- paste0(mutect_call," --dbsnp ", opt$dbSNP)
mutect_call <- paste0(mutect_call," --artifact_detection_mode ")


## on torquelord
} else { # Call variants genome - wide

				print("Parallelize variant calling")

				# Create tmp directory where to save temporary files for each chromosome
				dir.create(file.path(outvariants,caller,"tmp",sampleName),recursive = TRUE, showWarnings=FALSE)

				# create tmp directory on milton with all the chroms
				torque_tmp <- file.path("~",caller,"tmp")
				# concatenate and save it back

				chroms <- paste0("chr",c(1:22,"X","Y","M"))

				# create a script for every chromosome
				# I have the beamfiles in my unix directory
				# chrom script in my unix directory
				# save tmp and log files into torquelors mutect/tmp/

				for(chrom in chroms){

					mutect_call_chrom <- paste0(mutect_call," -L ",chrom,
						" -o ",file.path(torque_tmp,sampleName,paste0(type_call,"_snvs_indels_tmp_",chrom,".vcf")),
						" -log ", file.path(torque_tmp,sampleName,paste0(type_call,"_log_",chrom,".txt")))

					script <- file.path(outvariants,caller,"tmp",sampleName,paste0("script_",chrom,".sh"))
					cat("#!/bin/bash \n",file = script,append = FALSE)
					cat("module add gatk/3.7.0\n",file = script,append = TRUE)
					cat(mutect_call_chrom, file = script,append = TRUE)

					mutect_call_chrom <- NULL

				}

				##############################################################
				# Script to Concatenate variants called on separate chromosome
				##############################################################

				script_concat <- file.path(outvariants,caller,"tmp",sampleName,"concat_VCF.sh")
				concat_variants <- paste0("java -cp $GATK_JAR org.broadinstitute.gatk.tools.CatVariants", 
					" -R ",genome,
					paste0(" -V ",
						paste0(file.path(torque_tmp,sampleName),"/",
							paste0(type_call,"_snvs_indels_tmp_chr",c(1:2),".vcf"),collapse=" -V ")),
					" -out ",VCF," --log_to_file ",file.path(outvariants,caller,region_call,paste0(sampleName,"_log_catVariants")))

				cat("#!/bin/bash \n",file = script_concat,append = FALSE)
				cat("module add gatk/3.7.0\n",file = script_concat,append = TRUE)
				cat(concat_variants, file = script_concat,append = TRUE)

				#####################################
				# Create script to send to torquelord
				#####################################

				hpc_script <- file.path(outvariants,caller,"tmp",sample,"hpc_script.sh")

				cat("#!/bin/bash \n",file = hpc_script,append = FALSE)

				cat("ssh torquelord2.hpc.wehi.edu.au mkdir -p ",file.path(torque_tmp,sampleName),
					"\nfor chr in chr{1..2}"," ; do\nssh torquelord2.hpc.wehi.edu.au qsub -n -q small ",
					file.path(outvariants,caller,"tmp",sampleName,"script_${chr}.sh
					done\n"), 
					file = hpc_script,append = TRUE)

				system(paste0(file.path(outvariants,caller,"tmp",sample,"hpc_script.sh")))
				
				system(paste0("ssh torquelord2.hpc.wehi.edu.au qsub -n -q small ",script_concat))

				}

		}


module add gatk
gatk_nightly=/home/users/allstaff/quaglieri.a/software/GenomeAnalysisTK-nightly-2016-05-11-g3dbc89e/GenomeAnalysisTK.jar

parallel -j 2 gatk -T MuTect2 \
-R /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
--read_filter MappingQuality --min_mapping_quality_score 20 --max_alternate_alleles 4 \
-kmerSize 15 -kmerSize 30 \
-I:tumor /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/GRCh38_aligned/STAR_2pass/SRX322282_Recal.reorderedDupl.rg.split.bam \
--artifact_detection_mode  -L chr{}  \
-o /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/variant_calling/GRCh38/CD34/mutect/tmp/SRX322282_germline_snvs_indels_tmp_chr{}.vcf ::: {1..2}


-log /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/variant_calling/GRCh38/CD34/mutect/SRX322282_germline_log_chr{}.txt ::: {1..2}


###  Check why vep_parser does't work with the down_20M

local <- FALSE
dir <- ifelse(local,"/Volumes/AML_PROJECT/AML_RNA/cbf_aml_agrf","/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf")
mydir <- ifelse(local,"/Volumes/Anna's\ UNIX\ Home\ area/PHD_project/GEO_Leucegene_data","/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data")

# GATK and variant calling required files 
genome_fasta <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
genome_assembly <- "GRCh37"
bamfile <- "/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/hg19_aligned/STAR_2pass_20M_26880_AML/SRX729606Recal.reorderedDupl.rg.split.bam"
sample <- "SRX729606"
dbSNP <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/hg19_vcf/dbsnp_138.hg19.excluding_sites_after_129.vcf"
caller <- "mutect"
output_directory <- "/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/hg19_aligned/STAR_1pass_20M_26880_AML/../../variant_calling/hg19/down20M_26880"
regions <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/variant_calling/hg19/target_genes_Lavallee2016_unzip.bed"

# Rscript ${workdir}/call_variants.R --reference_fasta ${genome_fasta} --genome_assembly "GRCh37" \
# --bamfile ${bamin_final} --sampleName ${sample} --dbSNP ${dbSNP} -c "mutect" --output_directory ${outvariants} \
# --regions ${regions} 

opt <- data.frame(reference_fasta=genome_fasta,bamfile=bamfile,genome_assembly=genome_assembly,
	sampleName=sample,dbSNP=dbSNP,matched_normal=NA,
	variant_caller="mutect",output_directory=as.character(output_directory), create_PON=0, regions = regions)


18d17
< SRX729606
29d27
< SRX729619

# Normal sample
dir <- ifelse(local,"/Volumes/AML_PROJECT/AML_RNA/cbf_aml_agrf","/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf")
mydir <- ifelse(local,"/Volumes/Anna's\ UNIX\ Home\ area/PHD_project/GEO_Leucegene_data","/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data")

# GATK and variant calling required files 
genome_fasta <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
genome_assembly <- "GRCh37"
bamfile <- "/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/hg19_aligned/STAR_2pass_20M_26880_AML/SRX729629Recal.reorderedDupl.rg.split.bam"
sample <- "SRX729629"
dbSNP <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/hg19_vcf/dbsnp_138.hg19.excluding_sites_after_129.vcf"
caller <- "mutect"
output_directory <- "/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/hg19_aligned/STAR_1pass_20M_26880_AML/../../variant_calling/hg19/down20M_26880"
regions <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/variant_calling/hg19/target_genes_Lavallee2016_unzip.bed"

# Rscript ${workdir}/call_variants.R --reference_fasta ${genome_fasta} --genome_assembly "GRCh37" \
# --bamfile ${bamin_final} --sampleName ${sample} --dbSNP ${dbSNP} -c "mutect" --output_directory ${outvariants} \
# --regions ${regions} 

opt <- data.frame(reference_fasta=genome_fasta,bamfile=bamfile,genome_assembly=genome_assembly,
	sampleName=sample,dbSNP=dbSNP,matched_normal=NA,
	variant_caller="mutect",output_directory=as.character(output_directory), create_PON=0, regions = regions)


vep <- read.delim(annotatedVCF, fill = TRUE)
namesCol <- readLines(annotatedVCF,69)[[69]]
namesColsplit_w <- strsplit(as.character(namesCol),split="\t")[[1]]
colnames(vep_wrong) <- namesColsplit
colnames(vep_wrong)[1] <- "Uploaded_variation"

vep <- vep %>% tidyr::separate(Uploaded_variation, into = c("chrom","pos","genotype"), sep = "_")

print(paste0(basename(annotatedVCF_wrong), ": VEP call has been parsed."))

return(vep)

}


############
### VEP - VCF or TAB?
VCF <- "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/mutect/regions/66_germline_snvs_indels.vcf"
annotatedVCF <- "~/66_somatic_snvs_indels_annotated.vcf"

vep_call <- paste0("/usr/local/bioinfsoftware/VEP/VEP-85/bin/variant_effect_predictor.pl --dir_cache /usr/local/bioinfsoftware/VEP/VEP-85/.vep ",
				"-i ", VCF, " -o ",annotatedVCF, " --cache --everything --force_overwrite ", "--assembly GRCh38" ," --fork 12 --vcf")
system(vep_call)

bgzip /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/mutect/regions/6_somatic_snvs_indels.vcf
tabix -p vcf /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/mutect/regions/6_somatic_snvs_indels.vcf.gz
bgzip 6_somatic_snvs_indels_annotated.vcf
tabix -p vcf  6_somatic_snvs_indels_annotated.vcf.gz
vcf-merge ~/6_somatic_snvs_indels_annotated.vcf.gz \
/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/mutect/regions/6_somatic_snvs_indels.vcf.gz > ~/merged_6.vcf



echo -n "java -jar GenomeAnalysisTK.jar -T CombineVariants -R reference.fasta -o combined_output.vcf -genotypemergeOptions UNIQUIFY " && for vcf_file in $(ls folderWith100VcfFiles/*.vcf); do echo -n "--variant ${vcf_file}"; done"

#CHROM  POS     ID      REF     ALT     QUAL         FILTER                          
chr1    36466130        .       A       G       .       alt_allele_in_normal    

Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE"

AC=2;AN=8;CSQ=G|missense_variant|MODERATE|CSF3R|ENSG00000119535|Transcript|ENST00000331941|protein_coding|16/16||ENST00000331941.6:c.2318T>C|ENSP00000332180.5:p.Leu773Pro|2341|2318|773|L/P|cTt/cCt|||-1||SNV|HGNC|HGNC:2439||5||CCDS412.1|ENSP00000332180|Q99062||UPI000002AA5B|1|deleterious_low_confidence(0.01)|possibly_damaging(0.883)|||||||||||||||||||||||||||,G|3_prime_UTR_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000361632|protein_coding|15/15||ENST00000361632.8:c.*227T>C||2761|||||||-1||SNV|HGNC|HGNC:2439||1|P1|CCDS413.1|ENSP00000355406|Q99062||UPI000004CAC4|1|||||||||||||||||||||||||||||,G|3_prime_UTR_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000373103|protein_coding|17/17||ENST00000373103.5:c.*227T>C||3367|||||||-1||SNV|HGNC|HGNC:2439|YES|1||CCDS414.1|ENSP00000362195|Q99062||UPI000002AA5A|1|||||||||||||||||||||||||||||,G|missense_variant|MODERATE|CSF3R|ENSG00000119535|Transcript|ENST00000373104|protein_coding|18/18||ENST00000373104.5:c.2318T>C|ENSP00000362196.1:p.Leu773Pro|2866|2318|773|L/P|cTt/cCt|||-1||SNV|HGNC|HGNC:2439||1||CCDS412.1|ENSP00000362196|Q99062||UPI000002AA5B|1|deleterious_low_confidence(0.01)|possibly_damaging(0.883)|||||||||||||||||||||||||||,G|3_prime_UTR_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000373106|protein_coding|17/17||ENST00000373106.5:c.*227T>C||3286|||||||-1||SNV|HGNC|HGNC:2439||1|P1|CCDS413.1|ENSP00000362198|Q99062||UPI000004CAC4|1|||||||||||||||||||||||||||||,G|upstream_gene_variant|MODIFIER|MRPS15|ENSG00000116898|Transcript|ENST00000373116|protein_coding|||||||||||1693|-1||SNV|HGNC|HGNC:14504|YES|1|P1|CCDS411.1|ENSP00000362208|P82914||UPI0000135287||||||||||||||||||||||||||||||,G|upstream_gene_variant|MODIFIER|MRPS15|ENSG00000116898|Transcript|ENST00000462067|processed_transcript|||||||||||1785|-1||SNV|HGNC|HGNC:14504||2||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000464365|retained_intron|||||||||||2618|-1||SNV|HGNC|HGNC:2439||5|||||||1|||||||||||||||||||||||||||||,G|3_prime_UTR_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000464465|protein_coding|7/7||ENST00000464465.6:c.*227T>C||1393|||||||-1|cds_start_NF|SNV|HGNC|HGNC:2439||5|||ENSP00000435218||H0YE86|UPI0001F77A1D|1|||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000466138|retained_intron|||||||||||1099|-1||SNV|HGNC|HGNC:2439||3|||||||1|||||||||||||||||||||||||||||,G|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000480825|retained_intron|9/9||ENST00000480825.6:n.5988T>C||5988|||||||-1||SNV|HGNC|HGNC:2439||2|||||||1|||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000484762|retained_intron|||||||||||452|-1||SNV|HGNC|HGNC:2439||2|||||||1|||||||||||||||||||||||||||||,G|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000487540|processed_transcript|10/10||ENST00000487540.6:n.1919T>C||1919|||||||-1||SNV|HGNC|HGNC:2439||5|||||||1|||||||||||||||||||||||||||||;ECNT=1;HCNT=4;MAX_ED=.;MIN_ED=.;NLOD=287.72;SF=0f,1f;TLOD=7.30  GT:REF_F2R1:ALT_F1R2:QSS:AF:REF_F1R2:FOXOG:ALT_F2R1:AD  0/1:1:10:20309:0.014:645:0.00:0:646,10  0/0:1:29:44040:0.023:1410:0.00:0:1411,29        0/1:1:10:20309:0.014:645:0.00:0:646,10  0/0:1:29:44040:0.023:1410:0.00:0:1411,29

INFO    FORMAT  TUMOR   NORMAL  wehisan/general/user_managed/grpu_majewski_3/A
ML_RNA/cbf_aml_agrf/variant_calling/mutect/regions/6_somatic_snvs_indels_TUMOR  wehisan/general/user_managed/grpu_majewski_3/AML_RNA/c
bf_aml_agrf/variant_calling/mutect/regions/6_somatic_snvs_indels_NORMAL


######## VEP and VCF combine issues
# GATK and variant calling required files 

# Rscript ${workdir}/call_variants.R --reference_fasta ${genome_fasta} --genome_assembly "GRCh37" \
# --bamfile ${bamin_final} --sampleName ${sample} --dbSNP ${dbSNP} -c "mutect" --output_directory ${outvariants} \
# --regions ${regions} 

opt <- data.frame(reference_fasta=genome_fasta,bamfile=file.path(star_merged,cancer),
	genome_assembly="GRCh38",
	sampleName=sample,dbSNP=dbSNP,matched_normal=file.path(star_merged,normal),
	variant_caller="mutect",output_directory=outvariants,
	 create_PON=0, regions = regions)

annotatedVCF <- "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/varscan/regions/annotated_variants/30_somatic_annotated.vcf"


#############
### pdate VEP funciton in order to combine VEP and VCF in a sensible way

library(readr)
X6_somatic_snvs_indels_annotated <- read_delim("/Volumes/Anna's UNIX Home area/6_somatic_snvs_indels_annotated.vcf", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 238)
colnames(X6_somatic_snvs_indels_annotated) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO_VEP","FORMAT","TUMOR","NORMAL")
mu0_vcf$Location <- paste(mu0_vcf$chrom,mu0_vcf$pos,sep = "_")
X6_somatic_snvs_indels_annotated$Location <- paste(X6_somatic_snvs_indels_annotated$CHROM,X6_somatic_snvs_indels_annotated$POS,sep = "_")
merge_6 <- merge(mu0_vcf, X6_somatic_snvs_indels_annotated[,c("Location","INFO_VEP")])
dim(merge_6)

X66_somatic_snvs_indels_annotated <- read_delim("/Volumes/Anna's UNIX Home area/66_somatic_snvs_indels_annotated.vcf", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 236)
colnames(X66_somatic_snvs_indels_annotated) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO_VEP","FORMAT","SampleName")
mu0_vcf <- read.delim(file.path(dir,"variant_calling","mutect","regions/annotated_variants","66_germline_VCF_final.txt"),header = TRUE, fill = TRUE)
mu0_vcf$Location <- paste(mu0_vcf$chrom,mu0_vcf$pos,sep = "_")
X66_somatic_snvs_indels_annotated$Location <- paste(X66_somatic_snvs_indels_annotated$CHROM,X66_somatic_snvs_indels_annotated$POS,sep = "_")
merge_66 <- merge(mu0_vcf, X66_somatic_snvs_indels_annotated[,c("Location","INFO_VEP")])
dim(merge_66)

#X6_somatic_snvs_indels <- read_delim("/Volumes/AML_PROJECT/AML_RNA/cbf_aml_agrf/variant_calling/mutect/regions/6_somatic_snvs_indels.vcf", "\t", escape_double = FALSE, trim_ws = TRUE,skip = 236)
#colnames(X6_somatic_snvs_indels) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO_CALLER","FORMAT","TUMOR","NORMAL")
info_names <- as.character(read.delim("/Volumes/Anna's UNIX Home area/66_somatic_snvs_indels_annotated.vcf", skip = 235, nrow =1 , header = FALSE)[,1])
info_names1 <- gsub("##INFO=<ID=CSQ,Number=.,Type=String,Description=Consequence annotations from Ensembl VEP. Format: ","",info_names)
info_names1_string <- strsplit(info_names1, split = "|", fixed = TRUE)[[1]]


var <- merge_6[5,]

vep_expand <- function(var, VEP_Fields){
  
  info_vep_var <- var[names(var) == "INFO_VEP"]
  str1 <- gsub("=","",strsplit(as.character(info_vep_var), split = "CSQ", fixed = TRUE)[[1]][2])
  str2 <- strsplit(str1,split = ",", fixed = TRUE)[[1]]
  
  if ( length(str2) == 1 ){
    
    str1 <- paste0(str1,"|")
    str3 <- strsplit(str1,split = "|", fixed = TRUE)[[1]]
    names(str3) <- VEP_Fields
    
    info_variant <- var[1:(length(var) - 1)]
    info_variant_matrix <- data.frame(matrix(NA, nrow = 1, ncol = length(VEP_Fields) + length(info_variant)))
    colnames(info_variant_matrix) <- c(names(info_variant), VEP_Fields)
    info_variant_matrix[1,] <- c(info_variant, str3)
    final_vcf_vep <- info_variant_matrix
    
  } else {
    
    len <- length(str2)
    
    str3 <- lapply(str2,function(x) { 
      x <- paste0(x, "|")
      strsplit(x , split = "|", fixed = TRUE)[[1]]})
    
    str3_matrix <- do.call(rbind,str3)
    colnames(str3_matrix) <- VEP_Fields
    
    info_variant <- var[1:(length(var) - 1)]
    info_variant_matrix <- data.frame(matrix(info_variant, ncol = length(info_variant), nrow = len, byrow = TRUE))
    colnames(info_variant_matrix) <- names(info_variant)
    
    final_vcf_vep <- cbind(info_variant_matrix, str3_matrix)
    
    }
  
    return(vcf_vep = final_vcf_vep)

}


expand_variants <- do.call(rbind,apply(merge_6, 1 , vep_expand,  VEP_Fields = info_names1_string))

expand_variants$IMPACT_rank <- as.factor(expand_variants$IMPACT)
levels(expand_variants$IMPACT_rank) <- c(1,2,3,4)
expand_variants$IMPACT_rank <- as.character(as.numeric(expand_variants$IMPACT_rank))


multiqc ~/dir_bowtie/


## Debug luecegene data
# GATK and variant calling required files 
genome_fasta <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
genome_assembly <- "GRCh37"
bamfile <- "/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/hg19_aligned/STAR_2pass_30M_100_AML/SRX729607Recal.reorderedDupl.rg.split.bam"
sample <- "SRX729607"
dbSNP <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/hg19_vcf/dbsnp_138.hg19.excluding_sites_after_129.vcf"
caller <- "mutect"
output_directory <- "/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/hg19_aligned/STAR_1pass_20M_26880_AML/../../variant_calling/hg19/down30M_100"
regions <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/variant_calling/hg19/target_genes_Lavallee2016_unzip.bed"


opt <- data.frame(reference_fasta=genome_fasta,bamfile=bamfile,
	genome_assembly="GRCh37",
	sampleName=sample,dbSNP=dbSNP,matched_normal=NA,
	variant_caller="mutect",output_directory=output_directory,
	 create_PON=0, regions = regions)

#### STAR
star_genome=~/PHD_project/GEO_Leucegene_data/genomes/star_index_${assembly}
seed=100
reads="80M"
nreads=80000000
assembly="hg19"
workdir=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/quaglieri.a/scripts/functions
amldir=/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data
fastqdown_dir=${amldir}/sra_cbf_aml/fastq_${reads}_${seed}

Rscript ${workdir}/run_STAR.R --genome_index ${genome_fasta} \
--fastqfiles ${fastqdown_dir}/SRX729631_R1.fastq.gz,${fastqdown_dir}/SRX729631_R2.fastq.gz \
--sampleName SRX729631 --outdir ~/ 

#####################
###### VEP MAF vs AF
####################

vcf1=/Volumes/AML_PROJECT/quaglieri.a/GEO_Leucegene_data/variant_calling/hg19/down80M_100/varscan/regions/annotated_variants/SRX381851_germline_final.txt
vcf46=/Volumes/AML_PROJECT/quaglieri.a/GEO_Leucegene_data/variant_calling/hg19/down80M_100/varscan/regions/annotated_variants/SRX959064_germline_final.txt


#!/bin/bash

module load ensembl-vep 

set -e

vcf1=/Volumes/AML_PROJECT/quaglieri.a/GEO_Leucegene_data/variant_calling/hg19/down80M_100/varscan/regions/annotated_variants/SRX381851_germline_final.txt
vcf46=/Volumes/AML_PROJECT/quaglieri.a/GEO_Leucegene_data/variant_calling/hg19/down80M_100/varscan/regions/annotated_variants/SRX959064_germline_final.txt

vep --dir_cache /stornext/HPCScratch/cache/.vep -i $vcf1 \
-o ~/vcf1.vcf --cache --everything --force_overwrite --assembly GRCh38 --fork 12 --vcf

vep --dir_cache /stornext/HPCScratch/cache/.vep -i $vcf46 \
-o ~/vcf46.vcf --cache --everything --force_overwrite --assembly GRCh38 --fork 12 --vcf


qsub -n -q submit -l nodes=1:ppn=14,mem=20gb maf_af_vep.sh



#!/bin/bash


module load ensembl-vep

set -e

vcf1=/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/variant_calling/hg19/down80M_100/varscan/regions/SRX381851_germline_snvs_indels.vcf
vcf46=/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/variant_calling/hg19/down80M_100/varscan/regions/SRX959064_germline_snvs_indels.vcf

/usr/local/bioinfsoftware/VEP/VEP-85/bin/variant_effect_predictor.pl --dir_cache  /usr/local/bioinfsoftware/VEP/VEP-85/.vep -i $vcf1 \
-o ~/vcf1.vcf --cache --everything --force_overwrite --assembly GRCh38 --fork 12 --vcf

/usr/local/bioinfsoftware/VEP/VEP-85/bin/variant_effect_predictor.pl --dir_cache  /usr/local/bioinfsoftware/VEP/VEP-85/.vep -i $vcf46 \
-o ~/vcf46.vcf --cache --everything --force_overwrite --assembly GRCh38 --fork 12 --vcf

# on unix the columns are both GMAF ecc...

# torque new VEP version

#!/bin/bash

module load ensembl-vep

set -e

vcf1=/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/variant_calling/hg19/down80M_100/varscan/regions/SRX381851_germline_snvs_indels.vcf
vcf46=/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/variant_calling/hg19/down80M_100/varscan/regions/SRX959064_germline_snvs_indels.vcf

vep --dir_cache /stornext/HPCScratch/cache/.vep -i $vcf1 \
-o ~/vcf1.vcf --cache --everything --force_overwrite --assembly GRCh38 --fork 12 --vcf

vep --dir_cache /stornext/HPCScratch/cache/.vep -i $vcf46 \
-o ~/vcf46.vcf --cache --everything --force_overwrite --assembly GRCh38 --fork 12 --vcf

qsub -n -q submit -l nodes=1:ppn=14,mem=20gb maf_af_vep1.sh
qsub -n -q submit -l nodes=1:ppn=14,mem=20gb maf_af_vep46.sh


# now also torquelord gives the same output :/ but if writes AF instead of MAF :/:/:/ so I need to rerun
# all the others

# remove annotated variants removed for 20m_100


#############
### Exit code
#############

Rscript /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/quaglieri.a/scripts/functions/call_variants.R \
--reference_fasta /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
--genome_assembly "GRCh37" --bamfile /home/quaglieri.a/leucegene/down50_56745/STAR_2pass/SRX959064Recal.reorderedDupl.rg.split.bam \
--sampleName SRX959064 \
--dbSNP /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/hg19_vcf/dbsnp_138.hg19.excluding_sites_after_129.vcf \
-c mutect --output_directory /wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/hg19_aligned/STAR_1pass_50M_56745_AML/../../variant_calling/hg19/down50M_56745 \
--VEPcache /stornext/HPCScratch/cache/.vep --regions /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/variant_calling/hg19/target_genes_Lavallee2016_unzip.bed



VCF="/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/variant_calling/hg19/down30M_56745/mutect/regions/SRX729630_germline_snvs_indels.vcf"


#####
##### Whole genome variant calling debug
#####

script=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38_merged_runs/HPC_variant_calling/67_gatk_script_call_variants_varscan.sh
cp 67_gatk_script_call_variants_varscan.sh 67_gatk_script_call_variants_varscan_debug.sh


Rscript /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/../quaglieri.a/scripts/functions/call_variants_debug.R
 --reference_fasta /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg38/NCBI/Homo_sapiens/NCBI/GRCh38/Seq
uence/WholeGenomeFasta/genome.fa --bamfile /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh3
8_merged_runs/67_Recal.reorderedDupl.rg.split.bam --matched_normal /wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_
agrf/aligned_pass2_GRCh38_merged_runs/68_Recal.reorderedDupl.rg.split.bam --sampleName 67 -c varscan --output_directory /wehisan/
general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling

genome_fasta <- "/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
genome_assembly <- "GRCh37"
bamfile <- "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38_merged_runs/67_Recal.reorderedDupl.rg.split.bam"
matched_normal <- "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38_merged_runs/68_Recal.reorderedDupl.rg.split.bam"
sample <- "67"
caller <- "varscan"
output_directory <- "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling"

opt <- data.frame(reference_fasta=genome_fasta,
	bamfile=bamfile,
	genome_assembly="GRCh37",
	sampleName=sample,
	matched_normal=matched_normal,
	variant_caller=caller,
	output_directory=output_directory,
	create_PON=0,
	VEPcache="/usr/local/bioinfsoftware/VEP/VEP-85/.vep",
	Rrepos="https://cloud.r-project.org",
	RlibPath = "/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4")


### Debug STAR.R

# Get Counts with FeatureCounts 
local <- FALSE
dir <- ifelse(local,"/Volumes/AML_PROJECT/AML_RNA/cbf_aml_agrf","/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf")
mydir <- ifelse(local,"/Volumes/Anna's\ UNIX\ Home\ area/PHD_project/GEO_Leucegene_data","/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data")

# GATK and variant calling required files 
workdir <- file.path(dir, "../quaglieri.a/scripts/functions")
outdir <- file.path(dir,"aligned_pass2_GRCh38")
genome <- file.path("/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/star_index_GRCh38_99")
sampleName <- "91.ADE07KRH.PB.Dia_CAH6NANXX_L008"
STARmode <- "2PassMulti"
fastqdi <- "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/fastq/fastq_trimmed"
fastqfiles <- c(paste0(file.path(fastqdi,"8.ADE07KRH.PB.Rem_CAFMCANXX_TAATGCGC-TATAGCCT_L004_1P.fastq.gz"),",",
	file.path(fastqdi,"8.ADE07KRH.PB.Rem_CAFMCANXX_TAATGCGC-TATAGCCT_L004_2P.fastq.gz")))

bamout <- file.path(outdir,paste0(sampleName,"Aligned.sortedByCoord.out.bam"))

