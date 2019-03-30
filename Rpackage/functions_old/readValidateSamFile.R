#!/usr/bin/env Rscript

library(ggplot2)

## get arguments from command line
args <- commandArgs(trailingOnly=TRUE)

# validateDir <- "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38/ValidateSam"
validateDir <- args[1]
outdir <- args[2]

files <- list.files(validateDir,pattern="_validate.txt",full.names = TRUE)

readValidate <- function(validateSam=files){

	validateSam=file.path(validateSam)
	sampleName <- basename(validateSam)

	if(length(readLines(validateSam)) == 1){

		validateFile <- data.frame("No errors found",1,sampleName)
		colnames(validateFile) <- c("Error.Type","Count","sampleName")

	}else{
		validateFile <- data.frame(read.delim(validateSam,skip=3))
		validateFile$sampleName <- sampleName
	}

	return(validateFile)
}


errors <- do.call(rbind,lapply(files,readValidate))
colnames(errors)[1] <- "Error_type"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ggplot(data=errors,aes(x=sampleName,y=Count,fill=Error_type)) + geom_bar(position="stack",stat="identity")+coord_flip()
ggsave(file.path(outdir,"errorsBam_summary.pdf"),width=20,height=20)