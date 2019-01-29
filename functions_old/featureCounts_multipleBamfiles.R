#!/usr/bin/env Rscript

####################
# Read in arguments
####################

library(optparse)

print("Reading arguments in...")

option_list = list(
  
  make_option(c("--bamfiles"), type = "character", default = NULL, 
              help = "Path to file with a one column matrix containing path to bamfiles."),
 
  make_option(c("--annotation"), type = "character", default = "hg38", 
              help = "One of: hg19,hg38,mm10,mm9"),
  
  make_option(c("--runName"), type = "character", default = NULL, 
              help = "Run name"),
  
  make_option(c("--isPairedEnd"), type = "logical", default = TRUE, 
              help = "TRUE if bamfiles contains PE reads."),
  
  make_option(c("--Rrepos"), default = "https://cloud.r-project.org",
              help = "Redirection to server worldwide. Default: 'https://cloud.r-project.org' to install packages without setting a mirror."), 
  
  make_option(c("--RlibPath"), type = "character", default = "/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4",
              help = "R path to install R packages. Default: '/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4'.")
  
); 

# Parse arguments
opt_parser = OptionParser(option_list=option_list,add_help_option = TRUE);
opt = parse_args(opt_parser);

# Warning, error messages and sets defaults

# Check for Inconsistencies
if (!file.exists(opt$bamfiles)) {
  print_help(opt_parser)
  stop("No bamfile provided.", call.=FALSE)
}


#########################
# Load necessary packages
#########################


print("Loading R package and sourcing functions... ")

print(paste0("R repos : ", opt$Rrepos))

print(paste0("R lib path: ", opt$RlibPath))

repos_install <- opt$Rrepos

# Bioconductor packages
biocPackages <- c("Rsubread")
new.packages <- biocPackages[!(biocPackages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  source("https://bioconductor.org/biocLite.R")
  biocLite(new.packages,suppressUpdates = TRUE,suppressAutoUpdate = FALSE,ask = FALSE,siteRepos=opt$RlibPath)
}

lapply(biocPackages, require, character.only = TRUE,warn.conflicts = FALSE)


##############################
##### Get counts
#############################
# featureCounts
# CBF-AML

bamfiles <- read.table(opt$bamfiles,stringsAsFactors=FALSE)

# check esistance or run it
if(file.exists(file.path(dirname(bamfiles[1,1]),"featureCounts",
                         paste0("CBF-AML_hg19-geneCounts_Leucegene_FCinbuiltAnn_",opt$runName,".rds")))){
  
  print(paste0(
    file.path(dirname(bamfiles[1,1]),"featureCounts",
                         paste0("CBF-AML_hg19-geneCounts_Leucegene_FCinbuiltAnn_",opt$runName,".rds")),
    " already exists for this run"))

}else{

  # Using Inbuiltannotation

  fc <- featureCounts(files = bamfiles[,1],annot.inbuilt=opt$annotation, isPairedEnd = opt$isPairedEnd, nthreads = 15,ignoreDup=TRUE,minMQS=20)
  
  dir.create(file.path(dirname(bamfiles[1,1]),"featureCounts"),recursive = TRUE,showWarnings = FALSE)
  print(paste0("File saved in: ",file.path(dirname(bamfiles[1,1]),"featureCounts")))
  
  saveRDS(fc,file.path(dirname(bamfiles[1,1]),"featureCounts",paste0("CBF-AML_",opt$annotation,"-geneCounts_Leucegene_FCinbuiltAnn_",opt$runName,".rds")))

}


