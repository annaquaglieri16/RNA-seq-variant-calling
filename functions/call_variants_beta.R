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


####################
# Read in arguments
####################

print("Reading arguments in...")

# Useful when a matched normal is not available. In this way the variants called will be annotated with the PON.
# Define options for command line
option_list = list(
  
  make_option(c("-ref","--reference_fasta"), type = "character", default = NA, 
              help = "Path to the reference genome"),
  
  make_option(c("-assembly","--genome_assembly"), type = "character", default = "GRCh38", 
              help = "Genome assembly. One of: hg19 or GRCh38"),
  
  make_option(c("-bam", "--bamfile"), type = "character", default = NA, 
              help = "Full path to bamfiles"),
  
  make_option(c("-vcf", "--vcf"), type = "character", default = NA, 
              help = "Full path to vcf file."),
  
  make_option(c("-n", "--sampleName"), type = "character", default = NULL, 
              help = "Name for output files. Usually --bamfile without directory and extentions.", metavar = "character"),
  
  make_option(c("-c","--variant_caller"), type = "character", default = "mutect",
              help = "Variant caller to be used. Possible values are 'mutect', 'varscan', 'freebayes' and 'vardict'. Samtools mpileup is used to create VarScan input.",metavar = "character"),
  
  make_option(c("-type_call","--type_call"), type = "character", default = "germline",
              help = "germline or somatic call. Default is germline.",metavar = "character"),
  
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
              help = "If provided, limit the variant calling to the list of regions. Provide the list in the form of a bed file.",metavar = "character"),
  
  make_option(c("--create_PON"), type = "character", default = 0,
              help = "Possible vales: 0 and 1. Default 0. If create_PON = 1 option to create Panel of Normals are added."),
  
  make_option(c("--VEPcall"), type = "character", default = NA,
              help = "Call for VEP and cache directory. For example: vep --dir_cache /stornext/HPCScratch/cache/.vep"),
  
  make_option(c("--dbSNP"),type = "character",default = NULL,
              help = "dbSNP VCF for annotating variants with rs identifiers."),
  
  make_option(c("--Rrepos"), default = "https://cloud.r-project.org",
              help = "Redirection to server worldwide. Default: 'https://cloud.r-project.org' to install packages without setting a mirror."), 
  
  make_option(c("--RlibPath"), type = "character", default = "/stornext/HPCScratch/home/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.5",
              help = "R path to install R packages. Default: '/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4'."),
  
  make_option(c("--parse"), type = "logical", default = TRUE,
              help = ""),
  
  make_option(c("--force"), type = "logical", default = FALSE,
              help = "If TRUE then the variants are called, annotated and parsed even if they already exist.")
  
); 


# /usr/local/bioinfsoftware/VEP/VEP-85/.vep

# Parse arguments
opt_parser = OptionParser(option_list=option_list,add_help_option = TRUE);
opt = parse_args(opt_parser);

# Warning, error messages and sets defaults

# Check for Inconsistencies
if (!file.exists(opt$reference_fasta) & is.null(opt$vcf)) {
  print_help(opt_parser)
  stop("Reference genome is missing with no default.", call.=FALSE)
}

if (!file.exists(opt$bamfile) & is.null(opt$vcf)) {
  print_help(opt_parser)
  stop("Bamfile is missing with no default.", call.=FALSE)
}

#if (!file.exists(opt$bamfile) & opt$force) {
#  print_help(opt_parser)
#  stop("--force cannot be set to TRUE if bamfile is missing.", call.=FALSE)
#}

# Vardict cannot be run genome-wide

if (file.exists(opt$bamfile)) {
  if (is.null(opt$regions) & opt$variant_caller == "vardict") {
    stop("VarDict can only be run on a subset of the genome. Provide a bed file to the --regions argument.", call.=FALSE)
  }
}

# Warnings
# Use basename(bamfile) if the sample name is not provided
if (is.null(opt$sampleName)) {
  if(file.exists(opt$bamfile)){
    opt$sampleName <- gsub(".bam","",basename(opt$bamfile))
    warning("--sampleName is missing and basename(bamfile) will be used instead.", call. = TRUE)
    
  } else {
    
    opt$sampleName <- gsub(".vcf","",basename(opt$vcf))
    warning("--sampleName is missing and basename(vcf) will be used instead.", call. = TRUE)
  }
  
}

# Use the current directory as output directory unless otherwise stated
if (is.null(opt$output_directory)) {
  opt$output_directory <- getwd() 
  warning("--output_directory is missing. Files will be saved in the current directory.", call. = TRUE)
}


# No ref or bamfiles but there is a VCF file
if ((!file.exists(opt$reference_fasta) | !file.exists(opt$bamfile)) & !is.null(opt$vcf)) {
  warning("VCF file provided as input.", call.=TRUE)
}

#########################
# Load necessary packages
#########################


print("Loading R package and sourcing functions... ")

print(paste0("R repos : ", opt$Rrepos))

print(paste0("R lib path: ", opt$RlibPath))

repos_install <- opt$Rrepos

# R packages
list.of.packages <- c("tidyr","readr","doParallel","foreach","dplyr","devtools")
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

devtools::install_github("annaquaglieri16/samplepower")
library(samplepower)

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
vcf_input <- opt$vcf
call <- FALSE
annotate <- FALSE
parse <- FALSE

####################################
# Create prefix for output directory
####################################

# Regions or genome-wide

# Calling is performed only if a reference genome & bamfiles are provided as input
# If only VCF is given as input then the function will check if the annotated file exists.

# Define type of call: somatic or germline
type_call <- ifelse(file.exists(as.character(opt$matched_normal)), "somatic",opt$type_call)


if (file.exists(bamfile)) {
  
  region_call <- ifelse(!is.null(opt$regions) & is.na(opt$outputdir_suffix), "regions",
                        ifelse(!is.null(opt$regions) & !is.na(opt$outputdir_suffix), paste("regions",opt$outputdir_suffix,sep="_"),
                               ifelse(is.null(opt$regions) & is.na(opt$outputdir_suffix),"whole_genome",paste("whole_genome",opt$outputdir_suffix,sep="_"))))
  
  # Initialize log files
  # sink(file = file.path(outvariants,caller,region_call,paste0(sampleName,"_log_file_",Sys.Date())), type = "output")
  
}

print(Sys.time())

################
### Output files
################

# Create output directory
annot_directory <- ifelse(file.exists(vcf_input), 
                          file.path(dirname(vcf_input),"annotated_variants"),
                          file.path(outvariants,caller,region_call,"annotated_variants"))

# If there is a VCF file as input then save that one as VCF
VCF <- ifelse(file.exists(vcf_input),vcf_input,
              file.path(outvariants,caller,region_call,paste0(sampleName,"_",type_call,"_snvs_indels.vcf")))

logVCF <- ifelse(file.exists(vcf_input),gsub(".vcf","_log",vcf_input),
                 file.path(outvariants,caller,region_call,paste0(sampleName,"_",type_call,"_snvs_indels_log")))

annotatedVCF <- ifelse(file.exists(vcf_input),
                       file.path(dirname(vcf_input),"annotated_variants",paste0(sampleName,"_",type_call,"_annotated.vcf")),
                       file.path(outvariants,caller,region_call,"annotated_variants",paste0(sampleName,"_",type_call,"_annotated.vcf")))


parsedVCF <- ifelse(file.exists(vcf_input),
                    file.path(dirname(vcf_input),"annotated_variants",paste0(sampleName,"_",type_call,"_VCF_final.txt")),
                    file.path(outvariants,caller,region_call,"annotated_variants",paste0(sampleName,"_",type_call,"_VCF_final.txt")))


parsedVCF <- ifelse(file.exists(vcf_input),
                    file.path(dirname(vcf_input),"annotated_variants",paste0(sampleName,"_",type_call,"_VCF_final.txt")),
                    file.path(outvariants,caller,region_call,"annotated_variants",paste0(sampleName,"_",type_call,"_final.txt")))


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
  
  dir.create(annot_directory,recursive = TRUE,showWarnings = FALSE)
  
  
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
    
    dir.create(annot_directory,recursive = TRUE,showWarnings = FALSE)
    
  } else {
    
    # Check if parsedVCF exists
    options(warn=-1)
    check_existance_parsedVCF <- try(read.table(parsedVCF),silent=TRUE)
    options(warn=0)
      
    if ( class(check_existance_parsedVCF) == "try-error"  ) {
        
      print(paste0("Parse and extract fields from ", basename(annotatedVCF),"."))
        
      parse <- TRUE
        
    } else {
        
      print(paste0("VCF has already been annotated and parsed for sample ", sampleName,"."))
        
    }
  } 
}


# Check if annotation is needed anyway - for the moment a VCF file is parsed only if it gets annotated
annotate <- ifelse(is.na(opt$VEPcall), FALSE, annotate)
parse <- ifelse(is.na(opt$VEPcall), FALSE, parse)

# Force call, annotate and parsing
if(opt$force){
  
  call <- ifelse(file.exists(vcf_input), FALSE,TRUE)
  
  if(is.na(opt$VEPcall)){
    
    stop("No --VEPcall provided")
    
  }else{
    
    annotate <- TRUE
    parse <- TRUE
  }
}

print(paste0("Force: ", opt$force))
print(paste0("Caller: ", caller))
print(paste0("Call: ", call))
print(paste0("Annotate: ", annotate))
print(paste0("Parse: ", parse))

###############
# Call variants
###############

if (call) {
  
  print("Inside calling")
  
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
    
    
    if(is.null(opt$regions)){
      
      mutect_call <- paste0(mutect_call," -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 \
                            -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
                            -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 \
                            -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 \
                            -L chr21 -L chr22 -L chrX -L chrY -L chrM ", " -o ",VCF," -log ", logVCF)
      
    }
    
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


#parse <- ifelse(opt$parse,parse,opt$parse)

if ( parse ) {
  
  
  parse_vcf <- samplepower::parse_vcf_output(vcf_path = annotatedVCF,sample_name = sampleName, caller = caller)
  parse_vep_vcf <- samplepower::parse_vep_csq(vcf_path = annotatedVCF,vcf_df = parse_vcf)

  # Write final combined file to a txt file
  write.table(parse_vep_vcf,file = parsedVCF, quote = FALSE, sep = "\t",row.names = FALSE)
  
  print(paste0("Variant called and annotated for ", basename(VCF),"."))
  
}






