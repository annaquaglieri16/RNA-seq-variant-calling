#!/usr/bin/env Rscript

## To be tested!

# Modules needed: 
# module load STAR
# module load R

genome=$1
gtf=$2
genomeout=$3
genome_version=$4
maxSpan=$5


build_STAR_index <- function(genome_fasta,
							gtf,
							genomeout,
							genome_version = "hg38",
							maxSpan = 99){


# Checks
if (!file.exists(genome_fasta)) {
	stop("Reference genome is missing with no default.", call. = FALSE)
}

if (!file.exists(genome_fasta)) {
 	warning("GTF file is missing. Genome index will be built without information from gene annotation.",call. = TRUE)
}

if (!file.exists(genomeout)) {
	genomeout <- "~"
	warning("Genome index files will be output in the current directory.",call. = TRUE)
}


# Create output directory
genomedir <- file.path(genomeout,paste("star_index",genome_version,maxSpan,sep="_"))
dir.create(genomedir,recursive=TRUE,showWarnings=FALSE)


starIndex <- paste0("STAR --runThreadN $(nproc) --runMode genomeGenerate ",
"--genomeDir ${genomedir} \
--genomeFastaFiles ${genome} \
--sjdbOverhang ${maxSpan} \
--sjdbGTFfile ${gtf}")


print(starIndex)
system(starIndex)

}