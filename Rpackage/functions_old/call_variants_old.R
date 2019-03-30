#!/usr/bin/env Rscript

#########################
# Load necessary packages
#########################

library(optparse)

###############################################################
# Parse annotated VCF and allow to combine VarScan with Mutect2
###############################################################

# Take in input one file produced by vcf-query from VarScan and make it in a standard format
# Fields that I want
# Chrom pos ref alt vaf total_depth ref_depth alt_depth alt_for alt_rev ref_for ref_rev genotype genotype_quality caller

# geno_quality in varscan is ref_qual, alt_qual : phred score
# in mutect2 is the sum ref_qual+alt_qual
# in freebayes is ref_qual, alt_qual

################
## Varscan parser
################

# Get data
# Info fields from http://varscan.sourceforge.net/somatic-calling.html
# Simple example
# var1=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/variant_calling/varscan/regions/9_somatic_snvs_indels.snp.vcf
# vcf-query ${var1} -f '%CHROM %POS %REF %ALT %FILTER %INFO/SPV %INFO/GPV %INFO/SS [ %GT %GQ %DP %FREQ %RD %AD %DP4] \n' > ~/varscan_fields
# VCFvarscan <- "~/varscan_fields"

# only tumour
# # vcf-query ${var_old} -f '%CHROM %POS %REF %ALT %FILTER [ %GT %RBQ %ABQ %DP %FREQ %RD %AD %RDF %RDR %ADF %ADR] \n' > ~/var_old_fields
# VCFvarscan <- "~/var_old_fields"

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

			head(vcf)

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

			head(vcf_final)

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
	      vcf <- vcf %>% tidyr::separate(V9, into = c("ref_depth","alt_depth"))
	      vcf$ref_depth <- as.numeric(vcf$ref_depth)
	      vcf$alt_depth <- as.numeric(vcf$alt_depth)
	      vcf$tot_depth <- vcf$ref_depth + vcf$alt_depth
	      
	      # Allele frequency - and read alt_for/alt_rev - ref_for/ref_rev
	      vcf$V10 <- as.numeric(vcf$V10)
	      vcf$V11 <- as.numeric(vcf$V11)
	      vcf$V12 <- as.numeric(vcf$V12)
	      vcf$V13 <- as.numeric(vcf$V13)
	      vcf$V14 <- as.numeric(vcf$V14)
			
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

  # Freebayes only outputs AF as for a diploid organism. 
  # Therefore, in order to estimate a more refined VAF I need to take the alternative depth 
  # and the total depth. The AD field output from Freebayes provides the alternative depth information. 
  # However, it sometimes it outputs several values 2,3,4,5. By looking at IGV I decided to only 
  # consider the first of these values as well as the first value for alt_forw (SAF field) 
  # and alt_rev (SAR field). 

### parse freebayes
parse_freebayes_output <- function(VEPfreebayes){
  
  # Extract genotype
  readGeno <- readVcf(VEPfreebayes)
  genoDF <- data.frame(rownames(readGeno),
                       geno(readGeno)$GT,
                       geno(readGeno)$RO,
                       geno(readGeno)$AO)
  colnames(genoDF) <- c("MutName","genotype")
  
  
  
  # Convert VCF to a tsv table
  # module load vcflib
  system(paste0("vcf2tsv ",VEPfreebayes, " > " ,gsub(".vcf",".tsv",VEPfreebayes)))
  tsv.path <- gsub(".vcf",".tsv",VEPfreebayes)
  tsv <- read.delim(tsv.path)
  # Add genotype and caller columns
  tsv$MutName <- paste0(tsv$X.CHROM,":",tsv$POS,"_",tsv$REF,"/",tsv$ALT)
  tsv$Location <- paste(tsv$X.CHROM,tsv$POS,sep = "_")
  tsv <- merge(tsv,genoDF)
  tsv$caller <- "freebayes"
  tsv$MutName <- NULL
  
  # Subset columns to be consistent with the rest of the callers
  tsv_sub <- subset(tsv,select=c("Location","caller","X.CHROM","POS","REF","ALT","QUAL","FILTER","genotype",
                                 "DP","AF","RO","AO","QA","QR","SRF","SRR","SAF","SAR","CIGAR","CSQ"))
  
  
  tsv_sub <-  tsv_sub  %>% unite(geno_quality,QR,QA,sep=",")
  tsv_sub$alt_depth <- as.numeric(as.character(sapply(strsplit(as.character(tsv_sub$AO),split=","),function(allele) allele[[1]])))
  tsv_sub$alt_forw <- as.numeric(as.character(sapply(strsplit(as.character(tsv_sub$SAF),split=","),function(allele) allele[[1]])))
  tsv_sub$alt_rev <- as.numeric(as.character(sapply(strsplit(as.character(tsv_sub$SAR),split=","),function(allele) allele[[1]])))
  
  tsv_sub$VAF <- tsv_sub$alt_depth/(tsv_sub$alt_depth+tsv_sub$RO) # estimate of VAF
  tsv_sub$DP <- tsv_sub$alt_depth+tsv_sub$RO # estimate of total depth
  
  #plot(tsv_sub$alt_depth,tsv_sub$alt_forw+tsv_sub$alt_rev)
  #tsv_sub[tsv_sub$AO == "5,11,6",c("Location","X.CHROM","REF","ALT","AO","RO","DP","alt_depth","alt_forw","alt_rev")]
  #tsv_sub[tsv_sub$AO == "7,14,7,8",c("Location","X.CHROM","REF","ALT","AO","RO","DP","alt_depth","alt_forw","alt_rev")]
  
  tsv_sub$AF <- NULL
  tsv_sub$AO <- NULL
  tsv_sub$SAF <- NULL
  tsv_sub$SAR <- NULL
  colnames(tsv_sub) <- c("Location","caller","chrom","pos","ref","alt","qual","filter","genotype",
                         "tot_depth","ref_depth","geno_quality","ref_forw","ref_rev","cigar","INFO_VEP",
                         "alt_depth","alt_forw","alt_rev","VAF")
  
  tsv_final <- subset(tsv_sub,select = c("Location","caller","chrom","pos","ref","alt","qual","filter",
                                         "genotype","geno_quality","tot_depth","VAF","ref_depth","alt_depth","ref_forw","ref_rev","alt_forw","alt_rev","cigar","INFO_VEP"))
  
  return(tsv_final)
  
}



parse_vardict_output <- function(VEPvardict){
  
  # Extract genotype
  readGeno <- readVcf(VEPvardict)
  genoDF <- data.frame(rownames(readGeno),
                       geno(readGeno)$GT)
  colnames(genoDF) <- c("MutName","genotype")
  
  # Convert VCF to a tsv table
  # module load vcflib
  system(paste0("vcf2tsv ",VEPvardict, " > " ,gsub(".vcf",".tsv",VEPvardict)))
  tsv.path <- gsub(".vcf",".tsv",VEPvardict)
  tsv <- read.delim(tsv.path)
  # Add genotype and caller columns
  tsv$MutName <- paste0(tsv$X.CHROM,":",tsv$POS,"_",tsv$REF,"/",tsv$ALT)
  tsv$Location <- paste(tsv$X.CHROM,tsv$POS,sep = "_")
  tsv <- merge(tsv,genoDF)
  tsv$MutName <- NULL
  tsv$caller <- "vardict"
  
  # Subset columns to be consistent with the rest of the callers
  tsv_sub <- subset(tsv,select=c("Location","caller","X.CHROM","POS","REF","ALT","QUAL","FILTER","genotype",
                                 "DP","AF","ADJAF","VD","MQ","REFBIAS","VARBIAS","CSQ"))
  tsv_sub <-  tsv_sub  %>% 
    separate(REFBIAS,into=c("ref_forw","ref_rev"),sep=":") %>% 
    separate(VARBIAS,into=c("alt_forw","alt_rev"),sep=":") %>%
    mutate(ref_depth = DP - VD)
  
  colnames(tsv_sub) <- c("Location","caller","chrom","pos","ref","alt","qual","filter","genotype",
                         "tot_depth","VAF","ADJVAF_ADJ_indels","alt_depth","geno_quality",
                         "ref_forw","ref_rev","alt_forw","alt_rev","INFO_VEP","ref_depth")
  
  tsv_final <- subset(tsv_sub,select = c("Location","caller","chrom","pos","ref","alt","qual","filter",
                                         "genotype","geno_quality","tot_depth","VAF","ADJVAF_ADJ_indels","ref_depth","alt_depth","ref_forw","ref_rev","alt_forw","alt_rev","INFO_VEP"))
  
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
    if( length(grep("CSQ=",info_vep_var)) > 0 ){
      str1 <- gsub("=","",strsplit(as.character(info_vep_var), split = "CSQ=", fixed = TRUE)[[1]][2])
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

	make_option(c("-ref","--reference_fasta"), type = "character", default = NULL, 
		help = "Path to the reference genome"),

	make_option(c("-assembly","--genome_assembly"), type = "character", default = "GRCh38", 
		help = "Genome assembly. One of: hg19 or GRCh38"),

  	make_option(c("-bam", "--bamfile"), type = "character", default = NULL, 
              help = "Full path to bamfiles"),
	
	make_option(c("-n", "--sampleName"), type = "character", default = NULL, 
              help = "Name for output files. Usually --bamfile without directory and extentions.", metavar = "character"),
		
	make_option(c("-c","--variant_caller"), type = "character", default = "mutect",
		help = "Variant caller to be used. Possible values are 'mutect', 'varscan', 'freebayes' and 'vardict'. Samtools mpileup is used to create VarScan input.",metavar = "character"),

	make_option(c("-VD","--VarDict_dir"), type = "character", default = NA,
		help = "Path to VarDict extra functions.",metavar = "character"),

	make_option(c("--vaf"), type = "numeric", default = NA,
		help = "Variant allele frequency thesholds. For the moment only applied to VarDict.",metavar = "character"),

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

	make_option(c("--RlibPath"), type = "character", default = "/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.5",
		help = "R path to install R packages. Default: '/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4'."),

	make_option(c("--parse"), type = "logical", default = TRUE,
		help = "")

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

# Vardict cannot be run genome-wide
if (is.null(opt$regions) & opt$variant_caller == "vardict") {
  stop("VarDict can only be run on a subset of the genome. Provide a bed file to the --regions argument.", call.=FALSE)
}




#########################
# Load necessary packages
#########################


print("Loading R package and sourcing functions... ")

print(paste0("R repos : ", opt$Rrepos))

print(paste0("R lib path: ", opt$RlibPath))

repos_install <- opt$Rrepos

# R packages
list.of.packages <- c("tidyr","readr","doParallel","foreach","dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = repos_install, lib = opt$RlibPath)

lapply(list.of.packages, require, character.only = TRUE, warn.conflicts = FALSE,quietly=TRUE)

# Bioconductor packages
biocPackages <- c("VariantAnnotation")
new.packages <- biocPackages[!(biocPackages %in% installed.packages()[,"Package"])]

#source("https://bioconductor.org/biocLite.R")
#biocLite(suppressUpdates = TRUE,suppressAutoUpdate = FALSE,ask = FALSE)
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("IRanges","VariantAnnotation"),suppressUpdates = TRUE,suppressAutoUpdate = FALSE,ask = FALSE)

if(length(new.packages)){
	source("https://bioconductor.org/biocLite.R")
	biocLite(new.packages,suppressUpdates = TRUE,suppressAutoUpdate = FALSE,ask = FALSE,siteRepos=opt$RlibPath)
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

#####################################
# Define most used variables 
#####################################

genome <- opt$reference_fasta
bamfile <- opt$bamfile
caller <- opt$variant_caller
sampleName <- opt$sampleName
outvariants <- opt$output_directory
call <- FALSE
annotate <- FALSE
parse <- FALSE

####################################
# Create prefix for output directory
####################################

# Regions or genome-wide
region_call <- ifelse(!is.null(opt$regions) & is.na(opt$outputdir_suffix), "regions",
	ifelse(!is.null(opt$regions) & !is.na(opt$outputdir_suffix), paste("regions",opt$outputdir_suffix,sep="_"),
 ifelse(is.null(opt$regions) & is.na(opt$outputdir_suffix),"whole_genome",paste("whole_genome",opt$outputdir_suffix,sep="_"))))

# Define type of call: somatic or germline
type_call <- ifelse(file.exists(as.character(opt$matched_normal)), "somatic","germline")

# Create output directory
dir.create(file.path(outvariants,caller,region_call,"annotated_variants"),recursive = TRUE,showWarnings = FALSE)

# Initialize log files
# sink(file = file.path(outvariants,caller,region_call,paste0(sampleName,"_log_file_",Sys.Date())), type = "output")

print(Sys.time())

################
### Output files
################

VCF <- file.path(outvariants,caller,region_call,paste0(sampleName,"_",type_call,"_snvs_indels.vcf"))
logVCF <- file.path(outvariants,caller,region_call,paste0(sampleName,"_",type_call,"_snvs_indels_log"))
annotatedVCF <- file.path(outvariants,caller,region_call,"annotated_variants",paste0(sampleName,"_",type_call,"_annotated.vcf"))
fieldsVCF <- file.path(outvariants,caller,region_call,"annotated_variants",paste0(sampleName,"_",type_call,"_main_fields.txt"))
parsedVCF <- file.path(outvariants,caller,region_call,"annotated_variants",paste0(sampleName,"_",type_call,"_VCF_final.txt"))
combinedVepVcf <- file.path(outvariants,caller,region_call,"annotated_variants",paste0(sampleName,"_",type_call,"_final.txt"))

########################################
# CHECK EXISTANCE and decide step to run
########################################

# Check if VCF and/or annotate VCF exists
# If VCF does not exist: variant needs to be called (call=TRUE) and annotated (annotate = TRUE)
options(warn=-1)
check_existance_VCF <- try(read.table(VCF),silent=TRUE)
options(warn=0)

if(class(check_existance_VCF) == "try-error"){

	print(paste0("Call and annotate variants for sample ", sampleName)) 

	call <- TRUE
	annotate <- TRUE
	parse <- TRUE

} else {

	# If VCF exists: check existance of annotated VCF

	print(paste0("VCF file exists for sample ", sampleName))

	# Check if annotation exists
	options(warn=-1)
	check_existance_annotatedVCF <- try(read.table(annotatedVCF),silent=TRUE)
	options(warn=0)

	if (class(check_existance_annotatedVCF) == "try-error"){

		print(paste0("VCF files for sample ", sampleName, " needs to be annotated."))

		annotate <- TRUE 
		parse <- TRUE

	} else {

		if (caller == "mutect" | caller == "varscan"){

			# Check if annotation exists
			options(warn=-1)
			check_existance_fieldsVCF <- try(read.table(fieldsVCF),silent=TRUE)
			check_existance_parsedVCF <- try(read.table(parsedVCF),silent=TRUE)
			check_existance_combinedVepVcf <- try(read.table(combinedVepVcf),silent=TRUE)
			options(warn=0)

			if ( (class(check_existance_fieldsVCF) == "try-error") | (class(check_existance_combinedVepVcf) == "try-error") | (class(check_existance_parsedVCF) == "try-error")  ) {

				print(paste0("Parse and extract fields from ", basename(annotatedVCF),"."))

				parse <- TRUE

			} else {

				print(paste0("VCF has already been annotated and parsed for sample ", sampleName,"."))

			}

		} else { # for freebayes and vardict now
		
			options(warn=-1)
			check_existance_parsedVCF <- try(read.table(parsedVCF),silent=TRUE)
			check_existance_combinedVepVcf <- try(read.table(combinedVepVcf),silent=TRUE)
			options(warn=0)

			if ( (class(check_existance_parsedVCF) == "try-error") | (class(check_existance_combinedVepVcf) == "try-error") ) {

				print(paste0("Parse and extract fields from ", basename(annotatedVCF),"."))

				parse <- TRUE

			} else {

				print(paste0("VCF has already been annotated and parsed for sample ", sampleName,"."))

			} 
		}
	}
}

# Check if annotation is needed anyway
annotate <- ifelse(is.na(opt$VEPcall), FALSE, annotate)
parse <- ifelse(is.na(opt$VEPcall), FALSE, parse)

print(paste0("Caller: ", caller))
print(paste0("Call: ", call))
print(paste0("Annotate: ", annotate))
print(paste0("Parse: ", parse))

###############
# Call variants
###############

if (call) {

	if ( caller == "mutect" ) {

		mutect_call <- paste0("gatk -T MuTect2",
			" -R ", genome, 
			" -I:tumor ", bamfile)

		# Add matched normal if available
		if (file.exists(opt$matched_normal)){

			mutect_call <- paste0(mutect_call," -I:normal ", opt$matched_normal)

		}

		# Annotate with PON if available
		if (!is.null(opt$panel_of_normals)){
			# Cannot use it now cause it was done with hg19 and I need to redo it with GRCh38
			mutect_call <- paste0(mutect_call," --normal_panel ", opt$panel_of_normals)

		}

		# Add dbSNP annotation if available
		if (!is.null(opt$dbSNP)){

			mutect_call <- paste0(mutect_call," --dbsnp ", opt$dbSNP)

		}

		# Create Panel Of Normals?
		if ( opt$create_PON == 1){

			mutect_call <- paste0(mutect_call," --artifact_detection_mode ")

		}

		# Call on regions if a bed file is provided otherwise genome wide
		if (!is.null(opt$regions)){

			mutect_call <- paste0(mutect_call," -L ", opt$regions, " -o ",VCF," -log ", logVCF)

		} else { mutect_call <- paste0(mutect_call," -o ",VCF," -log ", logVCF) }

		print(mutect_call)
		system(mutect_call)
		print(paste0("Variant called with Mutect2 for ",basename(VCF),"."))

	}

	if ( caller == "varscan" ) {

		print(paste0("Calling variants for ",basename(VCF),"."))

		# Set preameter tp run VarScan based on type of call
		varscan_setting <- ifelse(type_call == "somatic","somatic","mpileup2cns")

		# Build samtools call
		samtools_call <- paste0("samtools mpileup ",
			"--output-tags AD,ADF,ADR,DP,SP ",
			"--fasta-ref ",genome)

		# Regions or genome wide
		if (!is.null(opt$regions)){

			samtools_call <- paste0(samtools_call," -l ",opt$regions," ")
		
		}

		# Decide if adding matched normal or tumour only
		if (type_call == "somatic"){ 

			# add normal tumour at the end of samtools_call
			samtools_call <- paste(samtools_call,opt$matched_normal,bamfile,sep=" ")
			
		} else { 

			# Tumour only call
			samtools_call <- paste(samtools_call,bamfile,sep=" ")

		}

		if (type_call == "somatic") {

			# Build varscan call
			varscan_call <- paste0("varscan ",varscan_setting, " -mpileup ", VCF,
			" --variants 1 --output-vcf 1 ")

		} else {

			# Build varscan call
			varscan_call <- paste0("varscan ",varscan_setting,
			" --variants 1 --output-vcf 1 --min-var-freq 0.01 > ",VCF) 

		}

			varscan_call <- paste0(varscan_call," 2> ",logVCF)

			# Combine samtools and varscan call
			samtools_varscan <- paste0(samtools_call," | ",varscan_call) 

			print(samtools_varscan)
			system(samtools_varscan)

			if (type_call == "somatic"){

				# When calling somatic calls Varscan produces two files .snp and .indel
				# ## Merge VCF so to have one VCF from the varscan output
				snp_out <- paste0(VCF, ".snp")
				indel_out <- paste0(VCF,".indel")
				
				system(paste("vcf-concat",snp_out,indel_out,">",VCF,sep=" "))

			}

			print(paste0("Variant called with VarScan for ",basename(VCF),"."))

	}

	########### FreeBayes
	## Only calls germline variants
	## Input should be bamfiles as out from STAR (sorted, markDupl, RG)

	if ( caller == "freebayes" ) {

	print(paste0("Calling variants for ",basename(VCF),"."))


	# Build samtools call
	freebayes_call <- paste0("freebayes --min-alternate-fraction 0.05 ",
		"--vcf ",VCF,
		" -f ",genome)

	# Regions or genome wide
	if (!is.null(opt$regions)){

		freebayes_call <- paste0(freebayes_call," --targets ",opt$regions," ")
	
	}

	# Only supports germline
	if (type_call == "somatic"){ 

		# add normal tumour at the end of samtools_call
		print("Freebayes only supports germline variant calling \n")
		freebayes_call <- paste0(freebayes_call," --bam ", bamfile)
		
	} else { 

		freebayes_call <- paste0(freebayes_call," --bam ", bamfile)

	}

		print(freebayes_call)
		system(freebayes_call)

		print(paste0("Variant called with freebayes for ",basename(VCF),"."))

	}

	## vardict can do somatic variants - filter duplicates by default
	# Default: 0x500 (filter 2nd alignments and duplicates).
	# performs by default indel realignment
	# -r 2 minimum n variant reads
	# Default coordinate system: 1 for BED file or amplicon BED file.

	if ( caller == "vardict" ) {

		print(paste0("Calling variants for ",basename(VCF),"."))

		vaf <- ifelse(is.na(opt$vaf),0.05,as.numeric(opt$vaf))
		# Build samtools call
		vardict_call <- paste0("vardict -f ", vaf, 
			" -c 1 -S 2 -E 3 -g 4 -r 2 -t -th 10 -v ",
			"-G ",genome)

		# Decide if adding matched normal or tumour only
		if (type_call == "somatic"){ 

			# add normal tumour at the end of samtools_call
			vardict_call <- paste0(vardict_call,
				" -b ",paste(bamfile,opt$matched_normal,sep="|"))
			
		} else { 

			# Tumour only call
			vardict_call <- paste(vardict_call,"-b",bamfile,sep=" ")

		}

		print(opt$VarDict_dir)
		# Regions are compulsory in VarDict 
		vardict_call <- paste0(vardict_call," ",opt$regions, " | ",
			file.path(opt$VarDict_dir,"teststrandbias.R"), " | ",
			file.path(opt$VarDict_dir,"var2vcf_valid.pl")," -N ",sampleName, " -E -f ",vaf," > " , VCF)

			print(vardict_call)
			system(vardict_call)

			print(paste0("Variant called with VarDict for ",basename(VCF),"."))

	}

} 



###################
# Annotate variants
###################

if (annotate) {

	#########
	# Run VEP 
	#########

	print(paste0("Annotate ",basename(VCF), " with VEP"))

	# output in tab format
	vep_call <- paste0(opt$VEPcall,
				" -i ", VCF, " -o ",annotatedVCF, " --cache --everything --force_overwrite ", "--assembly ", opt$genome_assembly ," --fork 12 --vcf")

	if (opt$genome_assembly == "GRCh37"){

		vep_call <- paste0(vep_call," --port 3337")

	}

	print(vep_call)
	system(vep_call)

}


parse <- ifelse(opt$parse,parse,opt$parse)

if ( parse ) {

	# to be changed for freebayes and vardict
	
	##############################
	# Extract fields from VCF call
	##############################

	if (caller == "mutect" | caller == "varscan") {

		print(paste0("Extract important fields from ",basename(VCF), " with vcf-query."))

		# VarScan looses three fields with tumour-only call
		fields_varscan <- ifelse(type_call == "somatic", "%INFO/SPV %INFO/GPV %INFO/SS [ %GT %GQ %DP %FREQ %RD %AD %DP4]", 
			"[ %GT %RBQ %ABQ %DP %FREQ %RD %AD %RDF %RDR %ADF %ADR]")

		# Mutect2 is the same for somatic and germline calls
		fields_mutect2 <- "[ %GT %QSS %AD %AF %ALT_F1R2 %ALT_F2R1 %REF_F1R2 %REF_F2R1]"

		# Final fields
		fields <- ifelse(caller == "mutect",fields_mutect2,fields_varscan)

		extract_fields <- paste0("vcf-query ", VCF, " -f '%CHROM %POS %REF %ALT %QUAL %FILTER ", fields," \n' > ", fieldsVCF)

		print(paste0("Extract fields from VCF from caller: ", caller))
		print(extract_fields)
		system(extract_fields)

		#############################################################################
		# Parsed extracted fields and create consistent output for VarScana and Mutect
		#############################################################################

		print(paste0("Parse ", basename(VCF), " and extract main fields."))

		if ( as.character(caller) == "varscan") {

			parsedVarscan <- parse_varscan_output(fieldsVCF, type_call = type_call)

			write.table(parsedVarscan,file = parsedVCF,quote = FALSE, sep = "\t",row.names = FALSE)

		} else {

			parsedMutect <- parse_mutect_output(fieldsVCF, type_call = type_call)

			write.table(parsedMutect,file = parsedVCF,quote = FALSE, sep = "\t",row.names = FALSE)

		}

	# Mutect and VaRscan parsed VCF files are the same
	}

	if ( caller == "freebayes"){

		parsedFreebayes <- parse_freebayes_output(annotatedVCF)

		write.table(parsedFreebayes,parsedVCF,col.names = TRUE,quote=FALSE,sep="\t")

	}

	if ( caller == "vardict"){

		parsedVardict <- parse_vardict_output(annotatedVCF)

		write.table(parsedVardict,parsedVCF,col.names = TRUE,quote=FALSE,sep="\t")

	}


	######################
	# Combine VEP and VCF
	######################

	# Using varianAnnotation package to extract exact fields
	get_info_fields_VEP <- readVcf(annotatedVCF)
	infos <- info(header(get_info_fields_VEP))
	info_names_string <- infos[rownames(infos) %in% "CSQ","Description"]
	info_names_string <- gsub("Consequence annotations from Ensembl VEP. Format: ","",info_names_string)
	info_names_string <- strsplit(info_names_string,split="[|]")[[1]]

	if (caller == "mutect" | caller == "varscan") {

		# Read parsed VCF from caller
		VCF_read <- read.delim(parsedVCF,header = TRUE, fill = TRUE)
		
		# Read VCF from VEP	
		skip_header <- (length(readLines(annotatedVCF)) - nrow(read.table(annotatedVCF)) )

		annotatedVCF_reads <- read_delim(annotatedVCF, "\t", escape_double = FALSE, trim_ws = TRUE, skip = skip_header - 1)
		#info_names <- as.character(read.delim(annotatedVCF, skip = skip_header - 2, nrow =1 , header = FALSE)[,1])
		#info_names1 <- gsub("##INFO=<ID=CSQ,Number=.,Type=String,Description=Consequence annotations from Ensembl VEP. Format: ","",info_names)
		#info_names_string <- strsplit(info_names1, split = "|", fixed = TRUE)[[1]]

		if (type_call == "somatic"){

			colnames(annotatedVCF_reads) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO_VEP","FORMAT","TUMOR","NORMAL") 
			
		} else {

			# colnames(annotatedVCF_reads)
	  		# [1] "#CHROM"  "POS"     "ID"      "REF"     "ALT"     "QUAL"    "FILTER"  "INFO"    "FORMAT"  "Sample1"
			colnames(annotatedVCF_reads) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO_VEP","FORMAT","TUMOR") 
		
		}

		# Add key identifier for the location of the variant
		VCF_read$Location <- paste(VCF_read$chrom,VCF_read$pos,sep = "_")
		annotatedVCF_reads$Location <- paste(annotatedVCF_reads$CHROM,annotatedVCF_reads$POS,sep = "_")

		# Merge VEP and VCF - still one line per variant
		merge_VCF_VEP <- merge(VCF_read, annotatedVCF_reads[,c("Location","INFO_VEP")], all.x = TRUE)
		
		} else { 

			# there is one field less now
			merge_VCF_VEP <- read.delim(parsedVCF)

	}

	# Combine VEP and VCF expanding variant information
	# This step is needed since sometimes the same variant called by a caller can have more than one 
	# annotation, say on different transcripts
	# I will create a data.frame where I repeat the same variant info from the caller (Chr, Pos, Allele) 
	# per annotation provided by VEP
	# expand_variants <- data.frame(do.call(rbind,apply(merge_VCF_VEP, 1 , vep_expand,  VEP_fields = info_names_string)))
	
	# Parallelize the process if more than 200 variants

	nvariants <- nrow(merge_VCF_VEP)
	maxcombine <- 10000

	if ( nvariants > 200 ) {

		no_cores <- detectCores() - 1
		# Initiate cluster
		cl <- makeCluster(no_cores)
		registerDoParallel(cl)

		seq_split <- seq(1,nvariants, by = 200)

		if ( seq_split[length(seq_split)] < nvariants ) seq_split <- c(seq_split,nvariants)

		#set.seed(100)
		#samplevar <- sample(nvariants,500)

		expand_variants <- foreach(nline = seq_split, 

		    .combine = rbind,
		    .inorder = FALSE, .maxcombine = maxcombine)  %dopar%  {

				start <- nline
				end <- ifelse(nline + 199 > nvariants, nvariants, nline + 199)

				data.frame(do.call(rbind,apply(merge_VCF_VEP[start:end,], 1, vep_expand ,
					VEP_fields = info_names_string)))

			}

		stopCluster(cl)
 
	} else {

		expand_variants <- data.frame(do.call(rbind,apply(merge_VCF_VEP, 1 , vep_expand,  VEP_fields = info_names_string)))

	}

	expand_variants$SampleName <- as.character(sampleName)

	# Make the IMPACT as number to be able to filter it out
	expand_variants$IMPACT <- as.character(expand_variants$IMPACT)
	expand_variants$IMPACT <- ifelse(expand_variants$IMPACT == "" , NA, expand_variants$IMPACT)
	expand_variants$IMPACT_rank <- as.factor(expand_variants$IMPACT)
	levels(expand_variants$IMPACT_rank) <- c(1,2,3,4)
	expand_variants$IMPACT_rank <- as.character(as.numeric(expand_variants$IMPACT_rank))


	# Write final combined file to a txt file
	write.table(expand_variants,file = combinedVepVcf, quote = FALSE, sep = "\t",row.names = FALSE)

	print(paste0("Variant called and annotated for ", basename(VCF),"."))

}



# sink()





