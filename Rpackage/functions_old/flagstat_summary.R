#!/usr/bin/env Rscript

#########################
# Load necessary packages
#########################

list.of.packages <- c("tidyr","optparse","dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)


####################
# Read in arguments
####################

option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="bamfiles directory", metavar="character"),
	make_option(c("-p", "--pattern"), type="character", default=NULL, 
              help="bamfile pattern to be searched for", metavar="character"),
	make_option(c("-o", "--outname"), type="character", default=NULL, 
              help="output name for combined csv file including all flagstats of this run.", metavar="character")
); 
 

# Parse arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$directory) | is.null(opt$pattern)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


# Assign command line arguments to R objects
dir <- opt$directory
pattern_file <- opt$pattern

dir <- "/wehisan/general/user_managed/grpu_majewski_3/quaglieri.a/GEO_Leucegene_data/hg19_aligned/STAR_2pass_AML"
dir2 <- "/wehisan/general/user_managed/grpu_majewski_3/TCGA_AML/bam_processed"
pattern_file <- "Aligned\\.reorderedDupl\\.rg\\.bam$"

data2  <- read.delim(file.path(dir2,"TCGA-AB-3012.PrimaryBloodDerivedCancer-PeripheralBlood.RNA.2b3.Dupl_flagstats"),header=FALSE)
data  <- read.delim(file.path(dir,"flagstats/SRX729607_flagstats"),header=FALSE)

# ?regex
bamfiles <- list.files(path = dir, pattern = pattern_file)

dir.create(file.path(dir,"flagstats"),recursive=TRUE,showWarnings=FALSE)

# Function to extract the properly mapped read pairs
getProperPairs <- function(bamfile){

	print(bamfile)
	sample_name <- gsub(pattern_file,"",bamfile)

	# If file is not readable then recompute sambamba flagstat
	outFile <- file.path(dir,"flagstats",
			paste0(sample_name,
				gsub(".bam","",bamfile),"_stats"))

	check_existance <- try(read.delim(outFile),
					 	silent=TRUE)

	if(class(check_existance) == "try-error"){

		# Compute flagstats for bamfile and save into a file
		system(paste0("sambamba flagstat ",file.path(dir,bamfile)," -t 5 > ",outFile))

		# Read the stats computed above and extract number of read pairs
		options(warn=-1)
		read_stats <- read.delim(outFile,stringsAsFactors=FALSE,header=FALSE)
		options(warn=0)
		colnames(read_stats) <- "stat"

		read_stats <- read_stats %>% 
		separate(stat,into=c("Count1","other"),sep="[+]",remove=FALSE) %>% 
		mutate(Count1 = gsub(" ","",Count1)) %>%
		separate(other,into=c("Count2","type"),sep=" ",remove=FALSE) 


		read_stats$Count1 <- gsub(" ","",read_stats$Count1)
		read_stats$Count2 <- gsub(" ","",read_stats$Count2)
		read_stats$Count1 <- as.numeric(as.character(read_stats$Count1))
		read_stats$Count2 <- as.numeric(as.character(read_stats$Count2))
		read_stats$TotReads <- read_stats$Count1 + read_stats$Count2

		prop_paired <- as.numeric(as.character(read_stats[grep("properly paired",read_stats$stat),"Count"]))

		return_count <- data.frame(ID=sample_name,properPaired=prop_paired)

		return(return_count) 

		} else {

		# If file already exists then read the stat file and extract the number of read pairs
		print(paste0("File ",outFile, " already exists." ))

		options(warn=-1)
		read_stats <- read.delim(outFile,stringsAsFactors=FALSE)
		options(warn=0)
		colnames(read_stats) <- "stat"
		read_stats <- read_stats %>% separate(stat,into=c("Count","other"),sep=" ",remove=FALSE)
		prop_paired <- as.numeric(as.character(read_stats[grep("properly paired",read_stats$stat),"Count"]))

		return_count <- data.frame(ID=sample_name,properPaired=prop_paired)

		return(return_count)

	}

}

# Get flagstat from all the bamfiles
proper_paired <- do.call(rbind,lapply(bamfiles,getProperPairs))

write.csv(proper_paired,file.path(dir,"flagstats",paste0("combine_flagstats_",opt$outname,".csv")),row.names=FALSE)

