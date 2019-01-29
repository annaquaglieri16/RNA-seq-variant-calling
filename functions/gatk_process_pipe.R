#!/usr/bin/env Rscript

# Modules needed: 
#module load sambamba
#module load gatk

library(optparse)

####################
# Read in arguments
####################

print("Reading arguments in...")

# Useful when a matched normal is not available. In this way the variants called will be annotated with the PON.
# Define options for command line
option_list = list(

	make_option(c("-ref","--reference_fasta"), type = "character", default = NULL, 
		help = "Path to the reference genome."),

  	make_option(c("-bam", "--bamfile"), type = "character", default = NULL, 
              help = "Full path to bamfile."),
	
	make_option(c("-n", "--sampleName"), type = "character", default = NULL, 
              help = "Name for output files. Usually --bamfile without directory and extentions.", metavar = "character"),
	
	make_option(c("--STARaligner"), type = "logical", default = TRUE, 
              help = "Logical: TRUE or FALSE, default to TRUE. If STAR has been used mapping qualities will be reassigned with GATK.", metavar = "character"),
	
	make_option(c("--knownSites1"),type = "character",default = NULL,
		help = "A database of known polymorphic sites. See GATK help page for BaseRecalibrator function."),

	make_option(c("--knownSites2"),type = "character",default = NULL,
		help = "A database of known polymorphic sites. See GATK help page for BaseRecalibrator function."),

	make_option(c("--knownSites3"),type = "character",default = NULL,
		help = "A database of known polymorphic sites. See GATK help page for BaseRecalibrator function."),

	make_option(c("--extraGATKoptions"),type = "character",default = NULL,
		help = "This options will be added to the GATK command when creating the Split.bam."),

	make_option(c("--RealignIndels"),type = "logical",default = FALSE,
		help = "Logical: TRUE or FALSE if GATK Indels realignment should be performed. Default to FALSE.")

#	make_option(c("--script"),type = "character",default = FALSE,
#		help = "Path where to save the script with the jobs to run."),

#	make_option(c("--run"),type = "logical",default = FALSE,
#		help = "Logical: TRUE or FALSE if the script with the jobs should be run.")

); 
 
# Parse arguments
opt_parser = OptionParser(option_list=option_list,add_help_option = TRUE);
opt = parse_args(opt_parser);

# Warning, error messages and sets defaults
# Error
if (!file.exists(opt$bamfile)) {
  print_help(opt_parser)
  stop("Bamfile is missing with no default.", call. = FALSE)
}

# Warnings
# Use basename(bamfile) if the sample name is not provided
if (is.null(opt$sampleName)) {
  opt$sampleName <- gsub(".bam","",basename(opt$bamfile))
  warning("--sampleName is missing and basename(bamfile) will be used instead.", call. = TRUE)
}


##########################################
# Provide info about the type of call made
##########################################

# List modules
print("Modules loaded")
system("module list")

####################################
# Create prefix for output directory
####################################

# Output directory - same as where the bamfile is saved
bamfile <- opt$bamfile
sampleName <- opt$sampleName
bamdir <- dirname(opt$bamfile)
suffix0 <- gsub(".bam","",basename(bamfile)) #bamfile without .bam
suffix <- gsub(sampleName,"",suffix0) # bamfile without .bam and SampleName

################
### Output files
################

# # File names
bamSplit <- file.path(bamdir,paste0(sampleName,suffix,".split.bam"))
bamSplitIndex <- file.path(bamdir,paste0(sampleName,suffix,".split.bam.bai"))

bamReal <- file.path(bamdir,paste0(sampleName,suffix,".split.realigned.bam"))
bamRealIndex <- file.path(bamdir,paste0(sampleName,suffix,".split.realigned.bam.bai"))

bamRecal <- file.path(bamdir,paste0(sampleName,suffix,".split.recalibrated.bam"))
bamRecalIndex <- file.path(bamdir,paste0(sampleName,suffix,".split.recalibrated.bam.bai"))

########################################
# CHECK EXISTANCE and decide step to run
########################################

check_existance_bamRecal <- try(file.info(bamRecal)$size, silent = TRUE)
check_existance_bamRecalIndex <- try(file.info(bamRecal)$size, silent = TRUE)
check_existance_bamRecal <- ifelse(is.na(check_existance_bamRecal),0,check_existance_bamRecal)
check_existance_bamRecalIndex <- ifelse(is.na(check_existance_bamRecalIndex),0,check_existance_bamRecalIndex)


splitNTrim <- FALSE
realignIndels <- FALSE
recalibrateBase <- FALSE
# scripts
splitCall <- NULL
splitBamIndexCall <- NULL
recalCall1 <- NULL
recalCall2 <- NULL
recalCall3 <- NULL
recalCall4 <- NULL
recalBamIndexCall <- NULL


if( (check_existance_bamRecal > 0) & (check_existance_bamRecalIndex > 0)){

	print(paste0(basename(bamRecal)," and ",basename(bamRecalIndex)," exist and are bigger than zero."))
	print(paste0("Call and annotate variants for sample ", sampleName))
	quit(save = "no") 

} else {

	# If bamRecal does not exist: check existance of bamReal
	check_existance_bamReal <- try(file.info(bamReal)$size, silent = TRUE)
	check_existance_bamRealIndex <- try(file.info(bamRealIndex)$size, silent = TRUE)
	check_existance_bamReal <- ifelse(is.na(check_existance_bamReal),0,check_existance_bamReal)
	check_existance_bamRecalIndex <- ifelse(is.na(check_existance_bamRealIndex) ,0,check_existance_bamRealIndex)


	if ( (check_existance_bamReal > 0) & (check_existance_bamRealIndex > 0)){

		print(paste0(basename(bamReal)," and ",basename(bamRealIndex)," exist and are bigger than zero."))

		recalibrateBase <- TRUE 

	} else {

		# If bamReal does not exist: check existance of bamSplit
		check_existance_bamSplit <- try(file.info(bamSplit)$size, silent = TRUE)
		check_existance_bamSplitIndex <- try(file.info(bamSplitIndex)$size, silent = TRUE)
		check_existance_bamSplit <- ifelse(is.na(check_existance_bamSplit),0,check_existance_bamSplit)
		check_existance_bamSplitIndex <- ifelse(is.na(check_existance_bamSplitIndex),0,check_existance_bamSplitIndex)

		if ( (check_existance_bamSplit > 0) & (check_existance_bamSplitIndex > 0) ) {

			print(paste0(basename(bamSplit)," and ",basename(bamSplitIndex)," exist and are bigger than zero."))

			realignIndels <- ifelse(realignIndels, TRUE, FALSE)
			recalibrateBase <- TRUE

		} else {

			splitNTrim <- TRUE
			realignIndels <- ifelse(opt$RealignIndels, TRUE, FALSE)
			recalibrateBase <- TRUE

		}
	}

}

print(paste0("Split'N'Trim: ", splitNTrim))
print(paste0("Realign Indels: ", realignIndels))
print(paste0("Recalibrate bases: ", recalibrateBase))

########################
### GATK pre-processing
########################

if (splitNTrim | realignIndels | recalibrateBase ){

	if (!file.exists(opt$reference_fasta)) {
  	print_help(opt_parser)
  	stop("Reference genome is missing with no default.", call.=FALSE)

	} else {

		genome_fasta <- opt$reference_fasta

	}

}

# Splin and Cigar
if(splitNTrim){

	splitCall <- paste0("gatk -T SplitNCigarReads ", 
		"-R ",genome_fasta,
		" -I ", bamfile,
		" -o ", bamSplit,
		" --filter_mismatching_base_and_quals -U ALLOW_N_CIGAR_READS")

	if ( opt$STARaligner){

		splitCall <- paste0(splitCall," -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 ")

	} 

	if (!is.null(opt$extraGATKoptions)){

		splitCall <- paste0(splitCall," ", opt$extraGATKoptions)

	} 

	splitCall <- paste0(splitCall," --log_to_file ", file.path(bamdir,paste0(sampleName,"_splitNtrim_log")))

	splitBamIndexCall <- paste0("sambamba index ",bamSplit)

}


if(realignIndels) {

	# to be written
	# realBamCall  

}


if(recalibrateBase){

	bamin <- ifelse(opt$RealignIndels, bamReal, bamSplit)

	# Create output directory
	baseRecalDir <- file.path(bamdir,"BaseQRecal",sampleName)
	dir.create(baseRecalDir, recursive = TRUE, showWarnings = FALSE)

	# Check presence of files for recalibration
	if (sum(file.exists(opt$knownSites1) & file.exists(opt$knownSites2) & file.exists(opt$knownSites3)) == 0 ) {
	  
	  print_help(opt_parser)
	  stop("At least one reference dataset should be provided to perform base reacalibration. See GATK BaseRecalibrator help page.", call.=FALSE)
	
	}

	knownSites <- c(opt$knownSites1,opt$knownSites2,opt$knownSites3)
	N_sites <- length(knownSites)
	knownSites_tmp <- paste0(rep(" -knownSites ", N_sites))
	knownSites_command <- paste(knownSites_tmp,knownSites, collapse = " ")

	# Step1 
	recalCall1 <- paste0("gatk -T BaseRecalibrator ",
		"-R ",genome_fasta,
		" -I ",bamin,
		" -nct 8",
		knownSites_command,
		" -o ",file.path(baseRecalDir,paste0(sampleName,"_recal_data.table")),
		" --log_to_file ", file.path(baseRecalDir,paste0(sampleName,"_recal_step1_log")))

	# Step2
	recalCall2 <- paste0("gatk -T BaseRecalibrator ",
		" -R ",genome_fasta,
		" -I ",bamin,
		" -nct 8",
		knownSites_command,
		" -BQSR ",file.path(baseRecalDir,paste0(sampleName,"_recal_data.table")),
		" -o ", file.path(baseRecalDir,paste0(sampleName,"_post_recal_data.table")),
		" --log_to_file ", file.path(baseRecalDir,paste0(sampleName,"_recal_step2_log")))

	# Step3
	recalCall3 <- paste0("gatk -T AnalyzeCovariates ",
		"-R ",genome_fasta,
		" -before ",file.path(baseRecalDir,paste0(sampleName,"_recal_data.table")),
		" -after ", file.path(baseRecalDir,paste0(sampleName,"_post_recal_data.table")),
		" -csv ", file.path(baseRecalDir,paste0(sampleName,"_recalibration_plots.csv")),
		" -plots ", file.path(baseRecalDir,paste0(sampleName,"_recalibration_plots.pdf")),
		" --log_to_file ", file.path(baseRecalDir,paste0(sampleName,"_recal_analyseCov_log")))

	# Step4 : recalibrate bases
	recalCall4 <- paste0("gatk -T PrintReads ",
		"-R ",genome_fasta,
		" -I ",bamin,
		" -o ", bamRecal,
		" -nct 8",
		" -BQSR ",file.path(baseRecalDir,paste0(sampleName,"_recal_data.table")),
		" --log_to_file ", file.path(baseRecalDir,paste0(sampleName,"_Log_recalibrated_bases")))

	recalBamIndexCall <- paste0("sambamba index ",bamRecal)

}


###############################################################
# Run jobs or only save scripts that then can be run on the HPC
###############################################################


print("Split and Trim GATK")
splitCall <- ifelse(splitNTrim,splitCall,"Split and trim already performed")
print(splitCall)
if(splitNTrim) system(splitCall)

print("")
print("Add index to Split and Trim bamfile")
splitBamIndexCall <- ifelse(splitNTrim,splitBamIndexCall,"Index created to split file")
print(splitBamIndexCall)
if(splitNTrim) system(splitBamIndexCall)


print("")
print("Recalibrate bases step 1")
recalCall1 <- ifelse(recalibrateBase,recalCall1,"Recal Step1 already performed")
print(recalCall1)
if(recalibrateBase) system(recalCall1)

print("")
print("Recalibrate bases step 2")
recalCall2 <- ifelse(recalibrateBase,recalCall2,"Recal Step2 already performed")
print(recalCall2)
if(recalibrateBase) system(recalCall2)

print("")
print("Recalibrate bases step 3")
recalCall3 <- ifelse(recalibrateBase,recalCall3,"Recal Step3 already performed")
print(recalCall3)
if(recalibrateBase) system(recalCall3)


print("")
print("Recalibrate bases step 4")
recalCall4 <- ifelse(recalibrateBase,recalCall4,"Recal Step4 already performed")
print(recalCall4)
if(recalibrateBase) system(recalCall4)

print("")
print("Add index to recalibrated bamfile")
recalBamIndexCall <- ifelse(recalibrateBase,recalBamIndexCall,"Recal Index already performed")
print(recalBamIndexCall)
if(recalibrateBase) system(recalBamIndexCall)


