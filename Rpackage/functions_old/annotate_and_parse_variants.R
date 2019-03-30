#!/usr/bin/env Rscript

#########################
# Load necessary packages
#########################

library(optparse)


getVAF_VarScan <- function(x){
	x <- as.character(x)
	x1 <- as.numeric(gsub("%","",x))/100
	return(x1)
}

separateDP4 <- function(x){

	x <- as.character(x)
	# Separate fields divided by a comma
	x1 <- 
	return(x1)
}


parse_varscan_output <- function(VCFvarscan, type_call = "somatic"){

	vcf <- read.table(VCFvarscan, fill = TRUE)

	if(type_call == "somatic"){

		if (dim(vcf)[2] < 22){

			stop(paste0(VCFvarscan, " has too few fields to come from a somatic call."))

		} else {

			print(paste0(VCFvarscan," are being extracted."))

			# 7 fields for Normal and 7 fields for Tumour

			vcf <- read.table(VCFvarscan)

			colnames(vcf) <- c("chrom","pos","ref","alt","qual","filter","somatic_pvalue","germline_pvalue","somatic_status",
				"genotype_N","geno_quality_N","tot_depth_N","freq_N","ref_depth_N","alt_depth_N","dp4_N",
				"genotype_T","geno_quality_T","tot_depth_T","freq_T","ref_depth_T","alt_depth_T","dp4_T")

			# Convert FREQ field in VarScan to a 0-1 VAF
			vcf$VAF_N <- do.call(c,lapply(vcf$freq_N,getVAF_VarScan))
			vcf$VAF_T <- do.call(c,lapply(vcf$freq_T,getVAF_VarScan))

			# Separate DP4 column
			vcf <- vcf %>% tidyr::separate(dp4_N, into = c("ref_forw_N","ref_rev_N","alt_forw_N","alt_rev_N"))
			vcf$ref_forw_N <- as.numeric(vcf$ref_forw_N)
			vcf$ref_rev_N <- as.numeric(vcf$ref_rev_N)
			vcf$alt_forw_N <- as.numeric(vcf$alt_forw_N)
			vcf$alt_rev_N <- as.numeric(vcf$alt_rev_N)
			vcf$ref_depth_N <- as.numeric(vcf$ref_depth_N)
			vcf$alt_depth_N <- as.numeric(vcf$alt_depth_N)

			vcf <- vcf %>% tidyr::separate(dp4_T, into = c("ref_forw_T","ref_rev_T","alt_forw_T","alt_rev_T"))
			vcf$ref_forw_T <- as.numeric(vcf$ref_forw_T)
			vcf$ref_rev_T <- as.numeric(vcf$ref_rev_T)
			vcf$alt_forw_T <- as.numeric(vcf$alt_forw_T)
			vcf$alt_rev_T <- as.numeric(vcf$alt_rev_T)
			vcf$ref_depth_T <- as.numeric(vcf$ref_depth_T)
			vcf$alt_depth_T <- as.numeric(vcf$alt_depth_T)

			# Select interestinf columns
			vcf$caller <- "varscan"

			vcf_final <- subset(vcf,select = c("caller","chrom","pos","ref","alt","qual","filter","somatic_pvalue","germline_pvalue","somatic_status",
				"genotype_N","geno_quality_N","tot_depth_N","VAF_N","ref_depth_N","alt_depth_N","ref_forw_N","ref_rev_N","alt_forw_N","alt_rev_N",
				"genotype_T","geno_quality_T","tot_depth_T","VAF_T","ref_depth_T","alt_depth_T","ref_forw_T","ref_rev_T","alt_forw_T","alt_rev_T"))

		}

	} else {

		if (dim(vcf)[2] > 17){

			stop(paste0(VCFvarscan, " has too many fields to come from a germline call."))

		} else {

			print(paste0(VCFvarscan," are being extracted."))

			vcf <- read.table(VCFvarscan)

			colnames(vcf) <- c("chrom","pos","ref","alt","qual","filter",
				"genotype","ref_base_quality","alt_base_quality","tot_depth",
				"freq","ref_depth","alt_depth","ref_forw","ref_rev","alt_forw","alt_rev")

			# Convert FREQ field in VarScan to a 0-1 VAF
			vcf$VAF <- do.call(c,lapply(vcf$freq,getVAF_VarScan))
			vcf$VAF <- as.numeric(vcf$VAF)
			
			#
			vcf$ref_forw <- as.numeric(vcf$ref_forw)
			vcf$ref_rev <- as.numeric(vcf$ref_rev)
			vcf$alt_forw <- as.numeric(vcf$alt_forw)
			vcf$alt_rev <- as.numeric(vcf$alt_rev)
			vcf$ref_depth <- as.numeric(vcf$ref_depth)
			vcf$alt_depth <- as.numeric(vcf$alt_depth)
			# Genotype quality definition - base quality from ref and rev
			vcf <- vcf %>% tidyr::unite(geno_quality,ref_base_quality,alt_base_quality,sep=",")

			# Select interestinf columns
			vcf$caller <- "varscan"

			vcf_final <- subset(vcf,select = c("caller","chrom","pos","ref","alt","qual","filter",
				"genotype","geno_quality","tot_depth","VAF","ref_depth","alt_depth","ref_forw","ref_rev","alt_forw","alt_rev"))

		}

	}

	return(vcf_final)

	print(paste0(VCFvarscan, ": Varscan VCF have been parsed."))

}

# parsed <- parse_varscan_output(VCFvarscan,"germline")

#################
## Mutect2 parser
#################

# Simple example
# mut=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/mutect/regions/30_somatic_snvs_indels.vcf
# vcf-query ${mut} -f '%CHROM %POS %REF %ALT %FILTER [ %GT %QSS %AD %AF %ALT_F1R2 %ALT_F2R1 %REF_F1R2 %REF_F2R1] \n' > ~/mutect_fields
# VCFmutect <- "~/mutect_fields"
# VCFmutect <- "~/mutect_old_fields"
# geno_quality_T = QSS = sum of base quality for each allele, probably QSS=ref,alt

parse_mutect_output <- function(VCFmutect, type_call = "somatic"){

	vcf <- read.table(VCFmutect, fill = TRUE)

	if(type_call == "somatic"){

		if (dim(vcf)[2] < 22){

			stop(paste0(VCFmutect, " has too few fields to come from a somatic call."))

		} else {

			print(paste0(VCFmutect," are being extracted."))

			# Ref and Alt depth 
			# mutect2 stores first tumour then normal
			vcf <- vcf %>% tidyr::separate(V8, into = c("ref_depth_T","alt_depth_T"))
			vcf$ref_depth_T <- as.numeric(vcf$ref_depth_T)
			vcf$alt_depth_T <- as.numeric(vcf$alt_depth_T)
			vcf$tot_depth_T <- vcf$ref_depth_T + vcf$alt_depth_T
			
			vcf <- vcf %>% tidyr::separate(V16, into = c("ref_depth_N","alt_depth_N"))
			vcf$ref_depth_N <- as.numeric(vcf$ref_depth_N)
			vcf$alt_depth_N <- as.numeric(vcf$alt_depth_N)
			vcf$tot_depth_N <- vcf$ref_depth_N + vcf$alt_depth_N

			# Allele frequency 
			# Tumour
			vcf$V9 <- as.numeric(vcf$V9)
			vcf$V10 <- as.numeric(vcf$V10)
			vcf$V11 <- as.numeric(vcf$V11)
			vcf$V12 <- as.numeric(vcf$V12)
			vcf$V13 <- as.numeric(vcf$V13)
			# Normal
			vcf$V17 <- as.numeric(vcf$V17)
			vcf$V18 <- as.numeric(vcf$V18)
			vcf$V19 <- as.numeric(vcf$V19)
			vcf$V20 <- as.numeric(vcf$V20)
			vcf$V21 <- as.numeric(vcf$V21)

			colnames(vcf) <- c("chrom","pos","ref","alt","qual","filter",
			"genotype_T","geno_quality_T","ref_depth_T","alt_depth_T", "VAF_T","alt_forw_T","alt_rev_T","ref_forw_T","ref_rev_T",
			"genotype_N","geno_quality_N","ref_depth_N","alt_depth_N","VAF_N","alt_forw_N","alt_rev_N","ref_forw_N","ref_rev_N","tot_depth_T","tot_depth_N")

			# Column to add to be consistent with VarScan
			vcf$caller <- "mutect2"
			vcf$somatic_pvalue <- NA
			vcf$germline_pvalue <- NA
			vcf$somatic_status <- vcf$filter

			vcf_final <- subset(vcf, select = c("caller","chrom","pos","ref","alt","qual","filter","somatic_pvalue","germline_pvalue","somatic_status",
				"genotype_N","geno_quality_N","tot_depth_N","VAF_N","ref_depth_N","alt_depth_N","ref_forw_N","ref_rev_N","alt_forw_N","alt_rev_N",
				"genotype_T","geno_quality_T","tot_depth_T","VAF_T","ref_depth_T","alt_depth_T","ref_forw_T","ref_rev_T","alt_forw_T","alt_rev_T"))

		}

	} else {

		# germline call 

		if (dim(vcf)[2] > 14){

			stop(paste0(VCFmutect, " has too many fields to come from a germline call."))

		} else {

			print(paste0(VCFmutect," are being extracted."))

			vcf <- read.table(VCFmutect)

			# Ref and Alt depth 
			vcf <- vcf %>% tidyr::separate(V8, into = c("ref_depth","alt_depth"))
			vcf$ref_depth <- as.numeric(vcf$ref_depth)
			vcf$alt_depth <- as.numeric(vcf$alt_depth)
			vcf$tot_depth <- vcf$ref_depth + vcf$alt_depth
			
			# Allele frequency 
			vcf$V9 <- as.numeric(vcf$V9)
			vcf$V10 <- as.numeric(vcf$V10)
			vcf$V11 <- as.numeric(vcf$V11)
			vcf$V12 <- as.numeric(vcf$V12)
			vcf$V13 <- as.numeric(vcf$V13)
			
			colnames(vcf) <- c("chrom","pos","ref","alt","qual","filter",
			"genotype","geno_quality","ref_depth","alt_depth","VAF","alt_forw","alt_rev","ref_forw","ref_rev","tot_depth")

			########
			# Column to add to be consistent with VarScan
			vcf$caller <- "mutect2"
			
			vcf_final <- subset(vcf,select = c("caller","chrom","pos","ref","alt","qual","filter",
				"genotype","geno_quality","tot_depth","VAF","ref_depth","alt_depth","ref_forw","ref_rev","alt_forw","alt_rev"))


		} 

		return(vcf_final)

		print(paste0(VCFmutect, ": mutect2 VCF have been parsed."))
	}

}



# # MUTECT call - mutect are the same instead VarScan change a bit : DP4 fields is 4 fields when germline call
# 	vcf-query ${mut} -f '%CHROM %POS %REF %ALT %FILTER [ %GT %QSS %AD %AF %ALT_F1R2 %ALT_F2R1 %REF_F1R2 %REF_F2R1] \n'

# 	vcf-query ${mut_old} -f '%CHROM %POS %REF %ALT %FILTER [ %GT %QSS %AD %AF %ALT_F1R2 %ALT_F2R1 %REF_F1R2 %REF_F2R1] \n'

# 	# VARSCAN calls
# 	vcf-query ${var1} -f '%CHROM %POS %REF %ALT %FILTER %INFO/SPV %INFO/GPV %INFO/SS [ %GT %GQ %DP %FREQ %RD %AD %DP4] \n'

# 	vcf-query ${var_old} -f '%CHROM %POS %REF %ALT %FILTER [ %GT %RBQ %ABQ %DP %FREQ %RD %AD %RDF %RDR %ADF %ADR] \n'


parse_freebayes_output <- function(VEPfreebayes){
  
  # Extract genotype
  readGeno <- readVcf(VEPfreebayes)
  genoDF <- data.frame(rownames(readGeno),
                       geno(readGeno)$GT)
  colnames(genoDF) <- c("MutName","genotype")
  
  # Convert VCF to a tsv table
  # module load vcflib
  system(paste0("vcf2tsv ",VEPfreebayes, " > " ,gsub(".vcf",".tsv",VEPfreebayes)))
  tsv.path <- gsub(".vcf",".tsv",VEPfreebayes)
  tsv <- read.delim(tsv.path)
  # Add genotype and caller columns
  tsv$MutName <- paste0(tsv$X.CHROM,":",tsv$POS,"_",tsv$REF,"/",tsv$ALT)
  tsv$MutName <- NULL
  tsv$Location <- paste(tsv$X.CHROM,tsv$pos,sep = "_")
  tsv <- merge(tsv,genoDF)
  tsv$caller <- "freebayes"
  
  # Subset columns to be consistent with the rest of the callers
  tsv_sub <- subset(tsv,select=c("Location","caller","X.CHROM","POS","REF","ALT","QUAL","FILTER","genotype",
                                 "DP","AF","RO","AO","QA","QR","SRF","SRR","SAF","SAR","CIGAR","CSQ"))
  tsv_sub <-  tsv_sub  %>% unite(geno_quality,QR,QA,sep=",")
  
  colnames(tsv_sub) <- c("Location","caller","chrom","pos","ref","alt","qual","filter","genotype",
                         "tot_depth","VAF","ref_depth","alt_depth","geno_quality","ref_forw","ref_rev","alt_forw","alt_rev","cigar","INFO_VEP")
  
  tsv_final <- subset(tsv_sub,select = c("Location","caller","chrom","pos","ref","alt","qual","filter",
                                         "genotype","geno_quality","tot_depth","VAF","ref_depth","alt_depth","ref_forw","ref_rev","alt_forw","alt_rev","INFO_VEP","cigar"))
 
  return(tsv_final)
  
}


################
## VEP parser
################

# Simple example
# /usr/local/bioinfsoftware/VEP/VEP-85/bin/variant_effect_predictor.pl --dir_cache /usr/local/bioinfsoftware/VEP/VEP-85/.vep \
# --assembly GRCh38 -i ${subset} -o ~/annotatedVCF_tab --cache --everything --force_overwrite --tab
# vcf-query ${vep_called} -f '%Uploaded_variation %Location %Allele %Gene %Feature_type %Consequence %Existing_variation [ %VARIANT_CLASS %IMPACT %SYMBOL ] \n' | head
# annotatedVCF="~/annotatedVCF_tab"

vep_parser <- function(annotatedVCF){

	vep <- read_delim(annotatedVCF, "\t", escape_double = FALSE, trim_ws = TRUE, skip = 68)
	colnames(vep)[1] <- "Uploaded_variation"

	vep <- vep %>% tidyr::separate(Uploaded_variation, into = c("chrom","pos","genotype"), sep = "_")

	print(paste0(basename(annotatedVCF), ": VEP call has been parsed."))

	return(vep)

}


##########################################
## Function to combine VEP with VCF output
##########################################

vep_expand <- function(var, VEP_fields){
  
 
  # var = one row = one variant from the VCF file
  # Extract INFO_VEP field from variant
  info_vep_var <- var[names(var) == "INFO_VEP"]
  # if no VEP annotation is present then there will be no |
  check_presence_variants <- try(grep("|",as.character(info_vep_var)))
  
  if( is.na(check_presence_variants) ){
    
    stop(paste0("VEP file ", info_vep_var," does not contain variants. Rerun VEP."), call.=FALSE)
    
  } else {
    
    # for mutect2 and varscan until I change them
    if( length(grep("CSQ",info_vep_var)) > 0 ){
      str1 <- gsub("=","",strsplit(as.character(info_vep_var), split = "CSQ", fixed = TRUE)[[1]][2])
    } else {
      str1 <- as.character(info_vep_var)
    }
    
	  # Chech if VEP output more than one annotation per variant
	  str2 <- strsplit(str1,split = ",", fixed = TRUE)[[1]]
	  
	  # If there is only one annotation from VEP than extract fields from INFO_VEP and concatenate them with the 
	  # Variant information
	  if ( length(str2) == 1 ){
	    
	    str1 <- paste0(str1,"|")
	    str3 <- strsplit(str1,split = "|", fixed = TRUE)[[1]]
	    names(str3) <- VEP_fields
	    
	    info_variant <- var[1:(length(var) - 1)]
	    info_variant_matrix <- data.frame(matrix(NA, nrow = 1, ncol = length(VEP_fields) + length(info_variant)))
	    colnames(info_variant_matrix) <- c(names(info_variant), VEP_fields)
	    info_variant_matrix[1,] <- c(info_variant, str3)
	    final_vcf_vep <- info_variant_matrix
	    
	  # If there is > 1 annotation per variant than repeat same information for the variant but with different annotations  
	  } else {
	    
	    len <- length(str2)
	    
	    str3 <- lapply(str2,function(x) { 
	      x <- paste0(x, "|")
	      strsplit(x , split = "|", fixed = TRUE)[[1]]})
	    
	    str3_matrix <- do.call(rbind,str3)
	    colnames(str3_matrix) <- VEP_fields
	    
	    info_variant <- var[1:(length(var) - 1)]
	    info_variant_matrix <- data.frame(matrix(info_variant, ncol = length(info_variant), nrow = len, byrow = TRUE))
	    colnames(info_variant_matrix) <- names(info_variant)
	    
	    final_vcf_vep <- cbind(info_variant_matrix, str3_matrix)
	    
    	}

    }
  
	# Return combined VCF and VEP   
    return(vcf_vep = final_vcf_vep)

}



####################
# Read in arguments
####################

print("Reading arguments in...")

# Useful when a matched normal is not available. In this way the variants called will be annotated with the PON.
# Define options for command line
option_list = list(


	make_option(c("-assembly","--genome_assembly"), type = "character", default = "GRCh38", 
		help = "Genome assembly. One of: hg19 or GRCh38"),
	
	make_option(c("-n", "--sampleName"), type = "character", default = NULL, 
              help = "Name for output files. Usually --bamfile without directory and extentions.", metavar = "character"),
		
	make_option(c("-c","--variant_caller"), type = "character", default = "mutect",
		help = "Variant caller to be used. Possible values are 'mutect', 'varscan', 'freebayes'. Samtools mpileup is used to create VarScan input.",metavar = "character"),

	make_option(c("--output_directory"), type = "character", default = NULL,
		help = "Output directory were VCF files are going to be saved."), 

	make_option(c("--outputdir_suffix"), type = "character", default = NA,
		help = "Suffix that will be postponed as mutect_outputdir_suffix."), 

	make_option(c("-m","--matched_normal"), type = "character", default = NA,
		help = "Path to the matched normal bamfile."),

	make_option(c("--panel_of_normals"), type = "character", default = NULL,
		help = "Path to a VCF file where variant have been called on a panel of normals."),

	make_option(c("--regions"), type = "character", default = NULL,
		help = "If provided, limit the variant calling to the list of regions provide the list in the form of a bamfile.",metavar = "character"),

	make_option(c("--create_PON"), type = "character", default = 0,
		help = "Possible vales: 0 and 1. Default 0. If create_PON = 1 option to create Panel of Normals are added."),

	make_option(c("--VEPcall"), type = "character", default = NA,
		help = "Call for VEP and cache directory. For example: vep --dir_cache /stornext/HPCScratch/cache/.vep"),

	make_option(c("--dbSNP"),type = "character",default = NULL,
		help = "dbSNP VCF for annotating variants with rs identifiers."),

	make_option(c("--Rrepos"), default = "https://cloud.r-project.org",
		help = "Redirection to server worldwide. Default: 'https://cloud.r-project.org' to install packages without setting a mirror."), 

	make_option(c("--RlibPath"), type = "character", default = "/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4",
		help = "R path to install R packages. Default: '/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4'.")

); 
 

# /usr/local/bioinfsoftware/VEP/VEP-85/.vep

# Parse arguments
opt_parser = OptionParser(option_list=option_list,add_help_option = TRUE);
opt = parse_args(opt_parser);

# Warning, error messages and sets defaults

# Check for Inconsistencies
if (!file.exists(opt$reference_fasta)) {
  print_help(opt_parser)
  stop("Reference genome is missing with no default.", call.=FALSE)
}

if (!file.exists(opt$bamfile)) {
  print_help(opt_parser)
  stop("Bamfile is missing with no default.", call.=FALSE)
}

# Warnings
# Use basename(bamfile) if the sample name is not provided
if (is.null(opt$sampleName)) {
  opt$sampleName <- gsub(".bam","",basename(opt$bamfile))
  warning("--sampleName is missing and basename(bamfile) will be used instead.", call. = TRUE)
}

# Use the current directory as output directory unless otherwise stated
if (is.null(opt$output_directory)) {
	opt$output_directory <- getwd() 
	warning("--output_directory is missing. Files will be saved in the current directory.", call. = TRUE)
}

#########################
# Load necessary packages
#########################


print("Loading R package and sourcing functions... ")

print(paste0("R repos : ", opt$Rrepos))

print(paste0("R lib path: ", opt$RlibPath))

repos_install <- opt$Rrepos

# R packages
list.of.packages <- c("tidyr","readr","doParallel","foreach")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = repos_install, lib = opt$RlibPath)

lapply(list.of.packages, require, character.only = TRUE, warn.conflicts = FALSE)

# Bioconductor packages
biocPackages <- c("VariantAnnotation")
new.packages <- biocPackages[!(biocPackages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
	source("https://bioconductor.org/biocLite.R")
	biocLite(new.packages,suppressUpdates = TRUE,suppressAutoUpdate = FALSE,ask = FALSE)
}

lapply(biocPackages, require, character.only = TRUE,warn.conflicts = FALSE)

####################################################
# Source function to parse output VCF and VEP call
####################################################

#source("/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/quaglieri.a/scripts/functions/parse_VCF_VEP_output.R")

print("Done!")

##########################################
# Provide info about the type of call made
##########################################

# print(paste0("Calling " , ifelse(!file.exists(opt$matched_normal), "somatic", "germline"), " mutations with ", opt$variant_caller," for sample ",opt$sampleName,"."))

# List modules
print("Modules loaded")
system("module list")
