#!/usr/bin/env Rscript

## To be tested!

# Modules needed: 
# module load STAR
# module load R

library(optparse)

####################
# Read in arguments
####################

print("Reading arguments in...")

# Useful when a matched normal is not available. In this way the variants called will be annotated with the PON.
# Define options for command line
option_list = list(

	make_option(c("--genome_index"), type = "character", default = NULL, 
		help = "Path to the folder with the reference genome index."),

  make_option(c("--fastqfiles"), type = "character", default = NULL, 
              help = "One or two comma separated full paths to the gzipped fastq files. \n If only one file is given STAR will consider it a SE library."),
	
	make_option(c("--sampleName"), type = "character", default = NULL, 
              help = "Name for output files. If not specified: --fastqfiles without directory and extention.", metavar = "character"),

	make_option(c("--outdir"), type = "character", default = NULL, 
              help = "Path to output directory. If not specified: ../STAR_align.", metavar = "character"),

	make_option(c("--sjfile"), type = "character", default = NULL, 
              help = "Path to output splije junction file from STAR 1-pass. Required if --STARmode '2PassMulti'.", metavar = "character"),
	
	make_option(c("--STARmode"), type = "character", default = "1Pass", 
              help = "One of: '2PassMulti', '2PassBasic', '1Pass'. For more information see the STAR manual for STAR 2-pass mode (Section 8)", metavar = "character"),

  make_option(c("--Rrepos"), default = "https://cloud.r-project.org",
    help = "Redirection to server worldwide. Need the default when installing packages without setting a mirror."), 

  make_option(c("--RlibPath"), type = "character", default = "/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4",
    help = "R path to install R packages.")

); 
 
# Parse arguments
opt_parser = OptionParser(option_list=option_list,add_help_option = TRUE);
opt = parse_args(opt_parser);

print(opt)

# Warning, error messages and sets defaults
# Error
if (!file.exists(opt$genome_index)) {
  print_help(opt_parser)
  stop("Reference genome is missing with no default.", call. = FALSE)
}

if (is.null(opt$fastqfiles)) {
  print_help(opt_parser)
  stop("FASTQ file is missing with no default.", call. = FALSE)
}

if ((opt$STARmode == "2PassMulti") & is.null(opt$sjfile)) {
  print_help(opt_parser)
  stop("--sjfile is required in --STARmode '2PassMulti'. See Section 8 of STAR Manual", call. = FALSE)
}

# Warnings
# Use basename(fastqfiles) if the sample name is not provided
if (is.null(opt$sampleName)) {
  opt$sampleName <- gsub(".fastq.gz","",basename(opt$fastqfiles[1]))
  warning("--sampleName is missing and basename(fastqfiles) will be used instead.", call. = TRUE)
}


# allow both PE and SE
# check if gzip otherwise gzip it


##########################################
# Provide info about the type of call made
##########################################

# print(paste0("Calling " , ifelse(!file.exists(opt$matched_normal), "somatic", "germline"), " mutations with ", opt$variant_caller," for sample ",opt$sampleName,"."))

# List modules
print("Modules loaded")
system("module list")


#########################
# Load necessary packages
#########################

print("Loading R package and sourcing functions... ")

print(paste0("R repos : ", opt$Rrepos))

print(paste0("R lib path: ", opt$RlibPath))

repos_install <- opt$Rrepos

# Bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsamtools")
require(Rsamtools)

list.of.packages <- c("parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = repos_install, lib = opt$RlibPath)

lapply(list.of.packages, require, character.only = TRUE)


#####################################
# Define most used variables 
#####################################

# Genome
genome <- opt$genome_index

# FQ files
fastqfiles <- opt$fastqfiles
fastqfiles1 <- gsub(" ","",strsplit(fastqfiles,split=",")[[1]])
nFQ <- length(fastqfiles1)
FQ1 <- fastqfiles1[1]
FQ2 <- ifelse(nFQ==2,fastqfiles1[2],NA)

if(nFQ > 2) warning("--fastqfiles should have at most two agruments. Onyl the first two will be used", call. = TRUE)

# Sample name
sampleName <- opt$sampleName

# Output directory
outdir <- opt$outdir

# STARmode
STARmode <- opt$STARmode

dir.create(opt$outdir,recursive=TRUE,showWarnings=FALSE)

# Print 
print(paste0("Genome :",genome))
print(paste0("FQ: ", FQ1, " ",FQ2))
print(paste0("Sample name: ",sampleName))
print(paste0("Output directory: ",outdir))
print(paste0("STAR mode: ",STARmode))

##################
# CHECK EXISTANCE 
##################

bamout <- file.path(outdir,paste0(sampleName,"Aligned.sortedByCoord.out.bam"))

options(warn=-1)
check_existance_BAM <- try(header <- scanBamHeader(bamout),silent=TRUE)
options(warn=0)

if(class(check_existance_BAM) == "try-error"){

  print(paste0("Align sample", sampleName)) 

  align <- TRUE

} else {

  print(paste0(basename(bamout)," has already been aligned."))

  align <- FALSE

}



############
# Run STAR 
############

if(align){

  # Check that the FQ files are gzipped otherwise compress them

  filetypeFQ1 <- summary( file(FQ1) )$class

  filetypeFQ2 <- ifelse(!is.na(FQ2),summary( file(FQ2) )$class, NA)

  if (filetypeFQ1 != "gzfile"){
    gzip(FQ1, destname=sprintf("%s.gz", FQ1), overwrite=FALSE, remove=TRUE)
    FQ1 <- sprintf("%s.gz", FQ1)
  }

  if(!is.na(filetypeFQ2) & filetypeFQ2 != "gzfile"){
    gzip(FQ2, destname=sprintf("%s.gz", FQ2), overwrite=FALSE, remove=TRUE)
    FQ2 <- sprintf("%s.gz", FQ2)
  }

  FasqtIn <- ifelse(!is.na(filetypeFQ2), paste(FQ1,FQ2),FQ1)

  # Build STAR call 

  no_cores <- detectCores() - 1

  STARcall <- paste0("STAR --genomeDir ",genome,
                     " --readFilesIn ", FasqtIn,
                     " --runThreadN ",no_cores,
                     " --chimSegmentMin 10 --readFilesCommand zcat --alignSJoverhangMin 8",
                     " --outBAMcompression 10 --alignSJDBoverhangMin 1 --limitBAMsortRAM 85741557872",
                     " --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 200000 --alignMatesGapMax 20000",
                     " --outFileNamePrefix ",file.path(outdir,sampleName),
                     " --outSAMtype BAM SortedByCoordinate")

# --outReadsUnmapped Fastx 

  if(STARmode == "1Pass"){

    STARcall <- paste0(STARcall, " --outFilterType BySJout --outFilterMultimapNmax 15")

  }


  if(STARmode == "2PassMulti"){

    STARcall <- paste0(STARcall, " --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3",
                        " --chimSegmentReadGapMax 6 --alignSJstitchMismatchNmax 5 -1 5 5 --chimOutType WithinBAM",
                        " --chimJunctionOverhangMin 2 --limitSjdbInsertNsj 2273673",
                        " --sjdbFileChrStartEnd ", as.character(opt$sjfile))

  }

  if(STARmode == "2PassBasic"){

    STARcall <- paste0(STARcall, " --twopassMode Basic --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3",
                        " --chimSegmentReadGapMax 6 --alignSJstitchMismatchNmax 5 -1 5 5 --chimOutType WithinBAM",
                        " --chimJunctionOverhangMin 2 --limitSjdbInsertNsj 2273673")

  }

  print(STARcall)
  system(STARcall)

}



