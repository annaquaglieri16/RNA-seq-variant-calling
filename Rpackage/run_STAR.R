#!/usr/bin/env Rscript

## To be tested!

# Modules needed: 
# module load STAR
# module load R

#library(optparse)

####################
# Read in arguments
####################

run_Star <- function(genome_index,
                      fastqfiles,
                      sampleName,
                      outdir,
                      sjfile,
                      STARmode = "2PassMulti",
                      Rrepos =  "https://cloud.r-project.org",
                      RlibPath = "/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4"){

# Warning, error messages and sets defaults
# Error
if (!file.exists(genome_index)) {
  stop("Reference genome is missing with no default.", call. = FALSE)
}

if (is.null(fastqfiles)) {
  stop("FASTQ file is missing with no default.", call. = FALSE)
}

if ((STARmode == "2PassMulti") & is.null(sjfile)) {
  stop("sjfile is required in STARmode '2PassMulti'. See Section 8 of STAR Manual", call. = FALSE)
}

# Warnings
# Use basename(fastqfiles) if the sample name is not provided
if (is.null(sampleName)) {
  sampleName <- gsub(".fastq.gz","",basename(fastqfiles[1]))
  warning("sampleName is missing and basename(fastqfiles) will be used instead.", call. = TRUE)
}


# allow both PE and SE
# check if gzip otherwise gzip it

print("Loading R package and sourcing functions... ")

print(paste0("R repos : ", Rrepos))

print(paste0("R lib path: ", RlibPath))

repos_install <- Rrepos

# Bioconductor
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rsamtools")
require(Rsamtools)

list.of.packages <- c("parallel","R.utils")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = repos_install, lib = opt$RlibPath)

lapply(list.of.packages, require, character.only = TRUE)

library(R.utils)
library(parallel)

#####################################
# Define most used variables 
#####################################

# Genome
genome <- genome_index

# FQ files
fastqfiles <- fastqfiles
fastqfiles1 <- gsub(" ","",strsplit(fastqfiles,split=","))
nFQ <- length(fastqfiles1)
FQ1 <- fastqfiles1[1]
FQ2 <- ifelse(nFQ==2,fastqfiles1[2],NA)

if(nFQ > 2) warning("fastqfiles should have at most two agruments. Onyl the first two will be used", call. = TRUE)

# Sample name
sampleName <- sampleName

# Output directory
outdir <- outdir

# STARmode
STARmode <- STARmode

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

  if(STARmode == "1Pass"){

    STARcall <- paste0(STARcall, " --outFilterType BySJout --outFilterMultimapNmax 15")

  }


  if(STARmode == "2PassMulti"){

    STARcall <- paste0(STARcall, " --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3",
                        " --chimSegmentReadGapMax 6 --alignSJstitchMismatchNmax 5 -1 5 5 --chimOutType WithinBAM",
                        " --chimJunctionOverhangMin 2 --limitSjdbInsertNsj 2273673",
                        " --sjdbFileChrStartEnd ", as.character(sjfile))

  }

  if(STARmode == "2PassBasic"){

    STARcall <- paste0(STARcall, " --twopassMode Basic --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3",
                        " --chimSegmentReadGapMax 6 --alignSJstitchMismatchNmax 5 -1 5 5 --chimOutType WithinBAM",
                        " --chimJunctionOverhangMin 2 --limitSjdbInsertNsj 2273673")

  }

  print(STARcall)
  system(STARcall)

}

}



