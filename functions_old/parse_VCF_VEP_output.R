###############################################################
# Parse annotated VCF and allow to combine VarScan with Mutect2
###############################################################

# Take in input one file produced by vcf-query from VarScan and make it in a standard format
# Fields that I want
# Chrom pos ref alt vaf total_depth ref_depth alt_depth alt_for alt_rev ref_for ref_rev genotype genotype_quality caller


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

		if (dim(vcf)[2] < 21){

			stop(paste0(VCFvarscan, " has too few fields to come from a somatic call."))

		} else {

			print(paste0(VCFvarscan," are being extracted."))

			# 7 fields for Normal and 7 fields for Tumour

			vcf <- read.table(VCFvarscan)

			colnames(vcf) <- c("chrom","pos","ref","alt","filter","somatic_pvalue","germline_pvalue","somatic_status",
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

			vcf_final <- subset(vcf,select = c("caller","chrom","pos","ref","alt","filter","somatic_pvalue","germline_pvalue","somatic_status",
				"genotype_N","geno_quality_N","tot_depth_N","VAF_N","ref_depth_N","alt_depth_N","ref_forw_N","ref_rev_N","alt_forw_N","alt_rev_N",
				"genotype_T","geno_quality_T","tot_depth_T","VAF_T","ref_depth_T","alt_depth_T","ref_forw_T","ref_rev_T","alt_forw_T","alt_rev_T"))

		}

	} else {

		if (dim(vcf)[2] > 16){

			stop(paste0(VCFvarscan, " has too many fields to come from a germline call."))

		} else {

			print(paste0(VCFvarscan," are being extracted."))

			vcf <- read.table(VCFvarscan)

			colnames(vcf) <- c("chrom","pos","ref","alt","filter",
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

			vcf_final <- subset(vcf,select = c("caller","chrom","pos","ref","alt","filter",
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

		if (dim(vcf)[2] < 21){

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

			colnames(vcf) <- c("chrom","pos","ref","alt","filter",
			"genotype_T","geno_quality_T","ref_depth_T","alt_depth_T", "VAF_T","alt_forw_T","alt_rev_T","ref_forw_T","ref_rev_T",
			"genotype_N","geno_quality_N","ref_depth_N","alt_depth_N","VAF_N","alt_forw_N","alt_rev_N","ref_forw_N","ref_rev_N","tot_depth_T","tot_depth_N")

			# Column to add to be consistent with VarScan
			vcf$caller <- "mutect2"
			vcf$somatic_pvalue <- NA
			vcf$germline_pvalue <- NA
			vcf$somatic_status <- vcf$filter

			vcf_final <- subset(vcf, select = c("caller","chrom","pos","ref","alt","filter","somatic_pvalue","germline_pvalue","somatic_status",
				"genotype_N","geno_quality_N","tot_depth_N","VAF_N","ref_depth_N","alt_depth_N","ref_forw_N","ref_rev_N","alt_forw_N","alt_rev_N",
				"genotype_T","geno_quality_T","tot_depth_T","VAF_T","ref_depth_T","alt_depth_T","ref_forw_T","ref_rev_T","alt_forw_T","alt_rev_T"))

		}

	} else {

		# germline call 

		if (dim(vcf)[2] > 13){

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
			
			colnames(vcf) <- c("chrom","pos","ref","alt","filter",
			"genotype","geno_quality","ref_depth","alt_depth","VAF","alt_forw","alt_rev","ref_forw","ref_rev","tot_depth")

			########
			# Column to add to be consistent with VarScan
			vcf$caller <- "mutect2"
			
			vcf_final <- subset(vcf,select = c("caller","chrom","pos","ref","alt","filter",
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
  # Remove initial fields which are not of interest
  check_presence_variants <- try(gsub("=","",strsplit(as.character(info_vep_var), split = "CSQ", fixed = TRUE)[[1]][2]))

  if( class(check_presence_variants) == "try-error" ){

  	stop(paste0("VEP file ", info_vep_var," does not contain variants. Rerun VEP."), call.=FALSE)

  } else {

	  str1 <- gsub("=","",strsplit(as.character(info_vep_var), split = "CSQ", fixed = TRUE)[[1]][2])
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

