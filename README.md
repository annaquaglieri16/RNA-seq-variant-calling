From FASQC to Variant Calling for RNA-Seq
================
Anna Quaglieri
11th October 2017

-   [Download RNA-Seq data from GEO](#download-rna-seq-data-from-geo)
    -   [Get SRX names](#get-srx-names)
    -   [Create NCBI query](#create-ncbi-query)
-   [Downsampling FASTQ files](#downsampling-fastq-files)
-   [Define files and programs needed for the pipeline](#define-files-and-programs-needed-for-the-pipeline)
-   [FASTQC](#fastqc)
-   [Alignement](#alignement)
    -   [Create `STAR` index](#create-star-index)
    -   [STAR-1pass](#star-1pass)
    -   [STAR-2pass](#star-2pass)
-   [GATK pre-processing](#gatk-pre-processing)
    -   [SplitNCigarReads](#splitncigarreads)
    -   [Base recalibration](#base-recalibration)
-   [Variant calling and annotation](#variant-calling)
-   [Used for `Finding a suitable library size for calling variants in RNA-Seq`](#used-for-finding-a-suitable-library-size-for-calling-variants-in-rna-seq)
    -   [Standardise the output of the annotated VCFs across callers](#standardise-the-output-of-the-annotated-vcfs-across-callers)
    -   [Depth of coverage for variants missed in downsampled runs](#depth-of-coverage-for-variants-missed-in-downsampled-runs)
    -   [Variant filtering](#variant-filtering)
        -   [Panel of normals (PON)](#panel-of-normals-pon)
        -   [Annotation databases and genomic features](#annotation-databases-and-genomic-features)
-   [Bibliography](#bibliography)

This is an example workflow from FASTQC to Variant calling using the functions in the *functions* folder.

``` bash
git clone git@github.com:annaquaglieri16/RNA-seq-variant-calling.git
```

All the functions used for the variant calling and downsampling pipeline are inside the `./functions` folder.

Download RNA-Seq data from GEO
------------------------------

### Get SRX names

``` r
library(GEOquery)
library(tidyverse)
library(knitr)
library(stringr)
```

Below is an example using one accession number from the Leucegene data.

``` r
dir.create("test_data",showWarnings = FALSE,recursive = TRUE)

# Get matrix files for every accession number
series_matrix_info <- function(gse){
gsed <- getGEO(gse,GSEMatrix=TRUE)
gse.mat <- pData(phenoData(gsed[[1]]))
reduced <- gse.mat[,c("title","geo_accession","relation.1")]
write.csv(reduced,file.path("test_data",paste(gse,"_",nrow(gse.mat),".csv",sep="")),row.names = FALSE)
}

series_matrix_info("GSE49642") # 43 samples
```

``` r
matrix_file <- list.files(path = file.path("test_data"),pattern = "GSE",full.names = TRUE)
GSEmatrix <- read_csv(matrix_file)

GSEmatrix$SRX <- str_extract(string = GSEmatrix$relation.1,pattern = "SRX[0-9][0-9][0-9][0-9][0-9][0-9]")
GSEmatrix$relation.1 <- NULL
kable(head(GSEmatrix))
```

| title  | geo\_accession | SRX       |
|:-------|:---------------|:----------|
| 02H053 | GSM1203305     | SRX332625 |
| 02H066 | GSM1203306     | SRX332626 |
| 03H041 | GSM1203307     | SRX332627 |
| 03H116 | GSM1203308     | SRX332628 |
| 03H119 | GSM1203309     | SRX332629 |
| 04H024 | GSM1203310     | SRX332630 |

### Create NCBI query

``` r
search_ncbi <- paste(GSEmatrix$SRX,collapse=" OR ")
search_ncbi
```

    ## [1] "SRX332625 OR SRX332626 OR SRX332627 OR SRX332628 OR SRX332629 OR SRX332630 OR SRX332631 OR SRX332632 OR SRX332633 OR SRX332634 OR SRX332635 OR SRX332636 OR SRX332637 OR SRX332638 OR SRX332639 OR SRX332640 OR SRX332641 OR SRX332642 OR SRX332643 OR SRX332644 OR SRX332645 OR SRX332646 OR SRX332647 OR SRX332648 OR SRX332649 OR SRX332650 OR SRX332651 OR SRX332652 OR SRX332653 OR SRX332654 OR SRX332655 OR SRX332656 OR SRX332657 OR SRX332658 OR SRX332659 OR SRX332660 OR SRX332661 OR SRX332662 OR SRX332663 OR SRX332664 OR SRX332665 OR SRX332666 OR SRX332667"

Paste the search SRX332625 OR SRX332626 OR SRX332627 OR SRX332628 OR SRX332629 OR SRX332630 OR SRX332631 OR SRX332632 OR SRX332633 OR SRX332634 OR SRX332635 OR SRX332636 OR SRX332637 OR SRX332638 OR SRX332639 OR SRX332640 OR SRX332641 OR SRX332642 OR SRX332643 OR SRX332644 OR SRX332645 OR SRX332646 OR SRX332647 OR SRX332648 OR SRX332649 OR SRX332650 OR SRX332651 OR SRX332652 OR SRX332653 OR SRX332654 OR SRX332655 OR SRX332656 OR SRX332657 OR SRX332658 OR SRX332659 OR SRX332660 OR SRX332661 OR SRX332662 OR SRX332663 OR SRX332664 OR SRX332665 OR SRX332666 OR SRX332667 into NCBI <https://www.ncbi.nlm.nih.gov/sra> and follow the intructions in <https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/#download-sequence-data-files-usi> **Download sequence data files using SRA Toolkit** to download all the `SRR` run names and information of the runs.

``` bash
# Files are saved in the home directory under ncbi/public/sra
prefetch --option-file SraAccList_CBF-AML_Leucegene.txt
```

SRA files were converted to `fastq` files with `fastq-dump --split-files`. The command `fastqc` was used to check the quality of the fastq files. The folder `test_data` contains small FASTQ files that can be used as example.

Downsampling FASTQ files
========================

The [seqtk](https://github.com/lh3/seqtk) tool was used to downsample an exact number of reads from paired end (PE) FASTQ files. The following is just an example command.

``` bash
seqtk sample -s100 test_data/SRR1608610_1.fastq.gz 10000 > test_data/sub_SRR1608610_1.fq
seqtk sample -s100 test_data/SRR1608610_2.fastq.gz 10000 > test_data/sub_SRR1608610_2.fq
```

Define files and programs needed for the pipeline
=================================================

-   The reference genome **hg19** is used for this analysis.
-   Below are all the programs and versions used

``` bash
module load STAR/2.5.2
module load R/3.4.3
module load anaconda2/4.0.0
module load sambamba/0.6.6
module load picard-tools/2.9.4
module load gatk/3.7.0
module load varscan/2.3.9
module load vcftools/0.1.13
module load samtools/1.6
module load ensembl-vep/89.0
module load vcflib/1.0.0-rc1
module load vardict/1.5.1
module load freebayes/1.1.0
module load picard-tools/2.9.4
```

-   The genome references and annotations used here have been downloaded from [iGenome website](https://support.illumina.com/sequencing/sequencing_software/igenome.html)

``` bash
# Hard link to genome.fa of the reference genome 
genome_fasta=path_to_hg19_genome_directory/genome.fa
# Hard link to gene.gtf where gene annotation is stored
gtf=path_to_hg19_gtf_directory/genes.gtf

# Functions directories
workdir=./functions

# folder where fastq files are stored
fastqdir=./test_data

# STAR folders for one-pass, two-pass and merged output
star_1pass=./aligned_star1
star_2pass=./aligned_star2
star_merged=./star_merged_runs # Every sample comes in different SRR runs which will have to be merged in one SRX sample.
```

FASTQC
======

This is just one example to run `fastqc` check on several FASTQ files.

``` bash
mkdir -p ${fastqdir}/fastqc
find ${fastqdir} -name "*.fastq.gz" > ./test_data/fastq_files.txt
cat ./test_data/fastq_files.txt | parallel -j 2 "fastqc {} --outdir ${fastqdir}/fastqc"
```

I strongly suggest to have a look at [MultiQC](http://multiqc.info/) which allows you to combine together the results of multiple samples into one compact document. To check which the programs whose output is supported by `MultiQC` and how to install it check their help page.

``` bash
multiqc ${fastqdir}/fastqc --interactive -n "FASTQC_summary" -o ${fastqdir}
```

The `FASTQC` reports offer a variety of measures and one can decide about discarding some samples or doing some adapter trimming if necessary. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) can be used for this purpose.

I suggest looking at one of my previous analysis around [adapters with STAR and Subread](https://github.com/annaquaglieri16/RNA-Seq-and-adapters--STAR-vs-Subjunc) since they can cause serious troubles with `STAR` default settings.

Alignement
==========

Once the *fastq* files are ready to be processed we can align them with *STAR*. [Subread/Rsubread](http://subread.sourceforge.net/) is another widely used RNA-Seq aligner. The reason why I initially choose *STAR* over *Subread* was simply due to the fact that *STAR* can generate specific output for chimeric reads that can be directly used with [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) to analyse gene fusions as well as because it is part of the [GATK Best Practices to call variants in RNA-Seq](https://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail).

Create `STAR` index
-------------------

`STAR` requires to build an index for the reference genome that will be used in the analysis. We are also going to detect fusions

``` bash
# TO BE DELETED
# Create Genome directory where to save STAR Index and STAR Fusion Index folders
star_genome100=path_to_genome_directory/star_index_hg19_99
star_fusion_data=path_to_genome_directory/star_fusion_hg19_dir

mkdir -p ${star_genome100}
mkdir -p ${star_fusion_data}
```

To build the *STAR index* one needs to provide the FASTA file for the reference genome used, a GTF file with information aabout the annotation and STAR also require an an extra parameter called `sjdbOverhang` which is usually set to be *(read length - 1)*. See `STAR` documentation for **Generating genome indexes** in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) - 99 is (read length - 1) as described in the first section.

Below is a wrapper for `STAR` call to build an index.

``` r
outpud_dir=./
${workdir}/build_STAR_index.sh ${genome_fasta} ${gtf} $outpud_dir "hg19" 99
```

STAR-1pass
----------

If you are working with a set of bamfiles, STAR developer suggests to run the alignment in a two-pass mode. This consists of first aligning all the bamfiles, collecting the *splice junctions* as output of STAR and realign all the bamfiles with this new information. For more details about 1-pass, 2-pass-multi and 2-pass-single see Section 8 of the [STAR documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). In my pipeline I normally use the 2-pass multi strategy as below.

``` bash
FQ1=test_data/SRR1608907_1.fastq.gz
FQ2=test_data/SRR1608907_2.fastq.gz

Rscript ./function/run_STAR.R --genome_index star_index_hg19_99 --fastqfiles $FQ1,$FQ2 \
    --sampleName SRR1608907 --outdir ./star_1pass --STARmode "1Pass" 
```

The `R` function above is a wrapper for the STAR call below:

``` bash
module add STAR/2.5

STAR --genomeDir path_to_star_index_hg19 \
--readFilesIn $FQ1 $FQ2 --runThreadN 27 --chimSegmentMin 10 --readFilesCommand zcat --alignSJoverhangMin 8 --outBAMcompression 10 --alignSJDBoverhangMin 1 --limitBAMsortRAM 85741557872 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 200000 --alignMatesGapMax 20000 --outFileNamePrefix path_to_star_1pass/SampleName --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax 15
```

To see all the arguments available:

``` bash
Rscript ./function/run_STAR.R --help
```

After running STAR on all the fastq files available we can collect all the splice junctions from the first pass and use them for the second pass.

``` bash
# concatenate splice junctions from all samples from ran in pass1
cat ./star_1pass/*SJ.out.tab > star_1pass/combined_sj.out.tab
# Dobin suggests to remove chrm cause they are usually False positives
awk '!/chrM/' ./star_1pass/combined_sj.out.tab > ./star_1pass/combined_sj_nochrM.out.tab
```

It is suggested to have a look at the summary of the alignment through *MultiQC*.

``` bash
multiqc star_1pass --interactive -n "STAR_1passQC" -o ./
```

STAR-2pass
----------

The second pass alignment is exactly the same as the first one with only a few differences:

-   the *sjfile* input created combining the splice junctions from the first pass
-   STAR is run with the option of output chimeric reads switched on. This will allow fusion analysis.

The ouput of STAR is a bamfile already sorted by coordinate with the suffix `Aligned.sortedByCoord.out.bam`. At this stage we can also run two more steps `post_align_qc1.sh` and `post_align_qc2.sh`.

``` bash
FQ1=test_data/SRR1608907_1.fastq.gz
FQ2=test_data/SRR1608907_2.fastq.gz

Rscript ./function/run_STAR.R --genome_index star_index_hg19_99 --fastqfiles $FQ1,$FQ2 \
    --sampleName SRR1608907 --outdir ./star_2pass --STARmode "2PassMulti" \
    --sjfile ./star_1pass/combined_sj_nochrM.out.tab

# Run featurecounts and collect fragment sizes for QC
./function/post_align_qc1.sh path_to_hg19_genome_directory/genome.fa \
path_to_hg19_gtf_directory/genes.gtf ./star_2pass/SRR1608907.Aligned.sortedByCoord.out.bam SRR1608907 

# Pre-process bamfile (add Read groups etc..)
./function/post_align_qc2.sh ./star_2pass/SRR1608907.Aligned.sortedByCoord.out.bam SRR1608907 \
path_to_hg19_genome_directory/genome.fa SRR1608907 >> ${script_dir}/${bamout}_align.sh
```

Details about the two post-alignment functions:

-   `post_align_qc1.sh` is optional:

1.  Runs [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) to get gene counts and compute PCA to evaulate the concordance between bamfiles sequenced on different lanes. This allows a QC before merging the different bamfiles into a single one.
2.  Runs [CollectMultipleMetrics](https://broadinstitute.github.io/picard/command-line-overview.html) to collect the fragment distribution of the bamfiles (only possible with PE reads). This is also a good QC to check that the fragment distribution of bamfiles on different lanes is the same.

-   `post_align_qc2.sh` contains necessary pre-prcessing steps:

1.  Marks PCR duplicates (using [sambamba markdup](http://lomereiter.github.io/sambamba/docs/sambamba-markdup.html))
2.  Add Read Groups to single runs before merging bamfiles (using [AddOrReplaceReadGroups](https://broadinstitute.github.io/picard/command-line-overview.html)). Even if files do not need to be merges, `GATK` requires read groups to be added bamfiles.
3.  Run [ValidateSamFile](https://broadinstitute.github.io/picard/command-line-overview.html) to check for errors in the final bamfile.

In order, its arguments are:

1.  **aligned bamfile**
2.  `SampleName`. This is the name of the sample applied to the `RGID` and `RGPU` fields below.
3.  `SampleName_run`. If a sample was sequenced across different lanes then these needs to be merged in one bamfile but lane-specific reag groups should be added to each separate file. This sample name will be used for the fields `RGLB` and `RGSM` in the `AddOrReplaceReadGroups` groups below.

``` bash
# Picard tool function to add read groups to a bamfile
AddOrReplaceReadGroups \
        I= ./star_2pass/SRR1608907.Aligned.sortedByCoord.out.bam \
        O= ./star_2pass/SRR1608907.Aligned.sortedByCoord.out.RG.bam \
        RGID=SRR1608907 \
        RGPU=SRR1608907 \
        RGLB=SRR1608907_run1 \
        RGPL="illumina" \
        RGSM=SRR1608907_run1
```

After running `post_align_qc2.sh` a file with the suffix `Aligned.reorderedDupl.rg.bam` will be created where read groups are added and PCR duplicated reads marked.

This time *MultiQC* will give us a summary output also of the fragment distributions and featureCounts if the output files are stored within the `star_2pass` folder.

``` bash
multiqc ./star_2pass --interactive -n "STAR_2passQC" -o ./
```

GATK pre-processing
===================

The pipeline contains function to call variants with `MuTect2`, `Samtools + VarScan2`, `VarDict` and `Freebayes`. In order to run `MuTect2` some GATK pre-processing are needed. The function `function/gatk_process_pipe.R` will perform the following steps:

-   *SplitNCigarReads* see [GATK documentation](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_rnaseq_SplitNCigarReads.php)
-   *Base recalibration* see [GATK documentation](https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr)

Below is an example call which wraps the steps above and check if files have already been created.

``` bash
Rscript ./functions/gatk_process_pipe.R \
--reference_fasta path_to_hg19_genome_directory/genome.fa \
--bamfile ./star_2pass/SRR1608907Aligned.reorderedDupl.rg.bam --sampleName SRX381851 \
-knownSites path_to_GATK_Bundle_files/dbsnp_138.hg19.excluding_sites_after_129.vcf \
-knownSites path_to_GATK_Bundle_files/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-knownSites path_to_GATK_Bundle_files/1000G_phase1.indels.hg19.sites.vcf 
```

The function above is a wrapper for the following `GATK3` calls.

SplitNCigarReads
----------------

``` bash
gatk -T SplitNCigarReads -R path_hg19_reference/genome.fa \
-I ./star_2pass/SRR1608907Aligned.reorderedDupl.rg.bam \
-o ./star_2pass/SRR1608907Aligned.reorderedDupl.rg.split.bam \
--filter_mismatching_base_and_quals -U ALLOW_N_CIGAR_READS -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
--log_to_file path_to_star_2pass/SRR1608907_RG_DUPL_SPLIT_log
```

Base recalibration
------------------

-   Base recalibration using known sites downloaded from the [GATK Bundle](https://github.com/snewhouse/ngs_nextflow/wiki/GATK-Bundle)

``` bash
module load gatk/3.7.0

gatk -T BaseRecalibrator -R path_hg19_reference/genome.fa \
-I ./star_2pass/SRR1608907Aligned.reorderedDupl.rg.split.bam -nct 8 \
-knownSites path_to_GATK_Bundle_files/dbsnp_138.hg19.excluding_sites_after_129.vcf \
-knownSites path_to_GATK_Bundle_files/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-knownSites path_to_GATK_Bundle_files/1000G_phase1.indels.hg19.sites.vcf \
-o ./star_2pass/BaseQRecal/SRR1608907/SRR1608907_recal_data.table \
--log_to_file ./star_2pass/BaseQRecal/SRR1608907/SRR1608907_recal_step1_log 

gatk -T BaseRecalibrator -R path_hg19_reference/genome.fa \
-I ./star_2pass/SRR1608907Aligned.reorderedDupl.rg.split.bam -nct 8 \
-knownSites path_to_GATK_Bundle_files/dbsnp_138.hg19.excluding_sites_after_129.vcf \
-knownSites path_to_GATK_Bundle_files/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-knownSites path_to_GATK_Bundle_files/1000G_phase1.indels.hg19.sites.vcf \
-BQSR path_to_star_2pass/BaseQRecal/SRR1608907/SampleName_recal_data.table \
-o path_to_star_2pass/BaseQRecal/SRR1608907/SRR1608907_post_recal_data.table \
--log_to_file path_to_star_2pass/BaseQRecal/SRR1608907/SRR1608907_recal_step2_log 

gatk -T AnalyzeCovariates -R path_hg19_reference/genome.fa \
-before path_to_star_2pass/BaseQRecal/SRR1608907/SSRR1608907_recal_data.table \
-after path_to_star_2pass/BaseQRecal/SRR1608907/SRR1608907_post_recal_data.table \
-csv path_to_star_2pass/BaseQRecal/SRR1608907/SSRR1608907_recalibration_plots.csv \
-plots path_to_star_2pass/BaseQRecal/SRR1608907/SRR1608907_recalibration_plots.pdf \
--log_to_file path_to_star_2pass/BaseQRecal/SRR1608907/SRR1608907_recal_analyseCov_log 

gatk -T PrintReads -R path_hg19_reference/genome.fa \
-I ./star_2pass/SRR1608907Aligned.reorderedDupl.rg.split.bam \
-o ./star_2pass/SRR1608907Recal.reorderedDupl.rg.split.bam \
-nct 8 -BQSR ./star_2pass/BaseQRecal/SRR1608907/SRR1608907_post_recal_data.table \
--log_to_file ./star_2pass/BaseQRecal/SRR1608907/SRR1608907_Log_recalibrated_bases
```

Variant calling and annotation
==============================

The target genes `./FindingSuitableCoverageData/Leucegene_target_genes.bed` were created using the gene symbols of the mutations listed in **Supplemental Table 3** of Lavallée VP et al 2016. We used the hg19 inbuilt annotation of `Rsubread` to obtain the gene ranges and added 500 bp at the end and at the beginning of each gene.

The function `call_variants.R` is a wrapper to call germline variants with MuTect2, VarScan2, VarDict and Freebayes and annotate them with `VEP`. The argument `--variant_caller` can be one of `mutect`, `varscan`, `freebayes` or `vardict`.

See all other options with:

``` bash
Rscript ./functions/call_variants.R
```

Below is an example to run MuTect2.

``` bash
Rscript ./functions/call_variants.R \
--reference_fasta path_hg19_reference/genome.fa \
--bamfile ./star_2pass/SRR1608907Recal.reorderedDupl.rg.split.bam --sampleName SRX959064 \
--regions ./FindingSuitableCoverageData/Leucegene_target_genes.bed --genome_assembly 'GRCh37' \
--variant_caller mutect --outputdir_suffix "caller_defaults" --VEPcall 'vep --dir_cache path_to_vep_cache_dir/.vep --offline' --output_directory ./variant_calling/
```

The call above is a wrapper for:

``` bash
gatk -T MuTect2 -R path_hg19_reference/genome.fa \
-I:tumor ./star_2pass/SRR1608907Recal.reorderedDupl.rg.split.bam \
-L ./FindingSuitableCoverageData/Leucegene_target_genes.bed \
-o ./variant_calling/mutect/region_caller_defaults/SRR1608907_germline_snvs_indels.vcf \
-log ./variant_calling/mutect/region_caller_defaults/SRR1608907_germline_snvs_indels_log
```

Or substituting `mutect` with `varscan` the following call will be made:

``` bash
samtools mpileup --output-tags AD,ADF,ADR,DP,SP \
--fasta-ref -R path_hg19_reference/genome.fa \
-l ./FindingSuitableCoverageData/Leucegene_target_genes.bed \
./star_2pass/SRR1608907Recal.reorderedDupl.rg.split.bam \
| varscan mpileup2cns --variants 1 --output-vcf 1 --min-var-freq 0.01 > \ ./variant_calling/varscan/region_caller_defaults/SRR1608907_germline_snvs_indels.vcf \ 
2> ./variant_calling/varscan/region_caller_defaults/SRR1608907_germline_snvs_indels_log
```

VarDict requires an extra argument which is a link to the path where the program was downloaded. This is because the program requires extra: `teststrandbias.R` and `var2vcf_valid.pl` scripts which were downloaded from <https://github.com/AstraZeneca-NGS/VarDict>.

``` bash
vardict -c 1 -S 2 -E 3 -g 4 -r 2 -t -th 10 -v -G \
-R path_hg19_reference/genome.fa \
-b ./star_2pass/SRR1608907Aligned.reorderedDupl.rg.bam \
./FindingSuitableCoverageData/Leucegene_target_genes.bed \
| vardict_dir/teststrandbias.R \
| vardict_dir/var2vcf_valid.pl -N -E -f 0.05 \
>  ./variant_calling/vardict/region_caller_defaults/SRR1608907_germline_snvs_indels.vcf
```

Below is an example using the output from VarScan but the same call is used for Mutect2 and VarDict vcf files.

``` bash
module load ensembl-vep/89.0

vep --dir_cache path_to_vep_cache_dir/.vep --offline \
-i ./variant_calling/varscan/region_caller_defaults/SRR1608907_germline_snvs_indels.vcf \
-o ./variant_calling/varscan/region_caller_defaults/annotated_variants/SRR1608907_germline_annotated.vcf \
--cache --everything --force_overwrite --assembly GRCh37 \
--fork 12 --vcf --port 3337
```

Used for `Finding a suitable library size for calling variants in RNA-Seq`
==========================================================================

Standardise the output of the annotated VCFs across callers
-----------------------------------------------------------

Since the aim of the paper `Finding a suitable coverage to call variants in RNA-Seq` is to compare the power in recovering variants with different library sizes and callers, we need first to standardize each caller's output. The VCF ouput from every caller is annotated with the Variant Effect Predictor version 89.0 and parsed within the function `call_variants.R` in order to produce a standardised and comparable output between callers. This output is not necessary is one only needs to analyse the VCF output from the calls above. Following, is the list of fields which are given as output after parsing the annotated VCF files. We provide a description of what information was extracted for every caller to populate each fields:

-   *Location*: `CHR_POS` where `CHR` and `POS` are standard VCF output of every caller.
-   *Caller*: either `mutect`, `varscan` or `vardict`.
-   *chrom*: field `CHROM` in the VCF file.
-   *pos*: field `POS` in the VCF file.
-   *ref*: field `REF` in the VCF file.
-   *alt*: field `ALT` in the VCF file.
-   *qual*: This field is not extracted consistently from the VCF output of the three callers but we made sure that its meaning should be consistent and it represents the average base quality at that position. Below is a description of how it is computed for each caller.
    -   In `MuTect2` the `QSS` field in the `FORMAT` fields reports the sum of the base qualities for the reference and alternative alleles separated by a comma. We used the reference and alternative depths at each position to compute the overall average base quality at that position. This quantity will populate the final `qual` field.
    -   In `VarScan2` this field is the average of the `RBQ` and `ABQ` fields from the `FORMAT` fields. They are defined respectively as the average quality of reference and alternative supporting bases in the header of the VCF file.
    -   In `VarDict` there are two `QUAL` fields, one is the standard 6th field of a VCF file and the other one is reported in the `INFO` fields. We will use the latter to populate the `qual` field in our analysis since it is the one defined as the average base quality at that position in the header of the VCF. VarDict uses a threshold of `QUAL >= 25` to report a variant.
-   *filter*: field `FILTER` in the VCF file. Each caller populates this field in different ways depending on the characteristics of the algorithm. In general, the entries for this field can be either `PASS` if that mutation passes all the filters determined by the caller or a description of the reason for filtering. The possible descriptions and their meaning can be found in the header of the VCF file generated by the caller.

-   *genotype*: field `GT` in the `FORMAT` fields of the VCF file.

-   *tot\_depth*: total read depth at each position as estimated by the caller. This information can be reported differently by each caller. `VarDict` and `VarScan` record it in the `DP` field while `MuTect2` records the reference and alternative depth in the `AD` columns and their sum was used to define the `total_depth`.

-   *VAF*: variant allele frequency for the variants recorded at that position. `VarScan` records it in the `FREQ` field while `MuTect` and `VarDict` in the `AF` field.

-   *ADJVAF\_ADJ\_indels*: field `ADJAF` only reported by `VarDict` and it represents the adjusted variant allele frequency for indels due to local realignment.

-   *ref\_depth*, *alt\_depth*, *ref\_forw*, *ref\_rev*, *alt\_forw* and *alt\_rev*: these fields represents the breakdown of supporting reference/alternative and forward/reverse reads at each location. Below is described what fields were used from every caller to extract these values. The fields are listed in order:
    -   `MuTect2`: *ref\_depth* and *alt\_depth* are the comma separated values reported in the field `AD`; *ref\_forw*, *ref\_rev*, *alt\_forw* and *alt\_rev* are respectively the `MuTect2` fields `REF_F1R2`, `REF_F2R1`, `ALT_F1R2` and `ALT_F2R1`.
    -   `VarScan2`: in order the features listed above are extracted from the the fields `RD`, `AD`, `RDF`, `RDR`, `ADF` and `ADR`.
    -   `VarDict`: *ref\_depth* and *alt\_depth* are the fields `REF` and `ALT` in the VCF file; the field `REFBIAS` contains comma separated values representing *ref\_forw* and *ref\_rev* and `VARBIAS` contains comma separated values representing *alt\_forw* and *alt\_rev*.

All the other fields are exactly as generated by the Variant Effect Predictor (ensembl-vep/89.0): `Allele`, `Consequence`, `IMPACT`, `SYMBOL`, `Gene`, `Feature_type`, `Feature`, `BIOTYPE`, `EXON`, `INTRON`, `HGVSc`, `HGVSp`, `cDNA_position`, `CDS_position`, `Protein_position` `Amino_acids`, `Codons`, `Existing_variation` `DISTANCE`, `STRAND`, `FLAGS`, `VARIANT_CLASS`, `SYMBOL_SOURCE`, `HGNC_ID`, `CANONICAL`, `TSL`, `APPRIS`, `CCDS`, `ENSP`, `SWISSPROT`, `TREMBL`, `UNIPARC`, `GENE_PHENO`, `SIFT`, `PolyPhen`, `DOMAINS`, `AF`, `AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF`, `SAS_AF`, `AA_AF`, `EA_AF`, `ExAC_AF`, `ExAC_Adj_AF`, `ExAC_AFR_AF`, `ExAC_AMR_AF`, `ExAC_EAS_AF`, `ExAC_FIN_AF`, `ExAC_NFE_AF`, `ExAC_OTH_AF`, `ExAC_SAS_AF`, `MAX_AF`, `MAX_AF_POPS`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `PUBMED`,`MOTIF_NAME`, `MOTIF_POS`, `HIGH_INF_POS`, `MOTIF_SCORE_CHANGE` `SampleName`, `IMPACT_rank`. Visit the VEP page to find more information <https://asia.ensembl.org/info/docs/tools/vep/index.html>.

Depth of coverage for variants missed in downsampled runs
---------------------------------------------------------

We used the GATK function `DepthOfCoverage` to recover the VAF, the total depth at a variant position and the depth of the alternative allele using the locations in (Lavallée et al. 2016).

``` bash
gatk -T DepthOfCoverage \
-R path_hg19_reference/genome.fa} \
-o ./SRR160890_loci_depth_of_coverage \
-I ./star_2pass/SRR1608907Recal.reorderedDupl.rg.split.bam \
--printBaseCounts \
-L ./FindingSuitableCoverageData/mutations_loci.bed 
```

Variant filtering
-----------------

Below is the code used to create the annotation information used to filter the variants called by the above pipeline (see *Methods* sections of <https://github.com/annaquaglieri16/Finding-a-suitable-library-size-for-calling-variants-in-RNA-Seq>).

### Panel of normals (PON)

17 CD34+ RNA-Seq samples (Accession Number GSE48846) were processed using the same pipeline reported above from download, to alignment and variant calling. It is important that the final VCF files used to create the Panel Of Normal (PON) reflect the pipeline used for the tumour samples. Therefore, the 14 normal libraries were aligned with STAR in the two pass mode and variants were called with VarDit, Samtools + VarScan2 and MuTect2 using the exact same code and pre-processing reported in the [Variant Calling section](#variant-calling). Each location `CHROM` and `POS` was used to filter variants found in the tumour samples. This implies that the PON variants didn't need to be annotated with VEP 89.0.

A PON was created for every caller using the code below. We are assuming that the variants have been called for all the CD34+ RNA-Seq samples and stored in `normal_cd34_caller`.

``` bash
module load gatk/3.7.0
module load vcflib/1.0.0-rc1

gatk -T CombineVariants \
-R path_hg19_reference/genome.fa \
--variant:SRX322282 normal_cd34_caller/SRX322282_germline_snvs_indels.vcf \
--variant:SRX322283 normal_cd34_caller/SRX322283_germline_snvs_indels.vcf \
--variant:SRX322284 normal_cd34_caller/SRX322284_germline_snvs_indels.vcf \
--variant:SRX322285 normal_cd34_caller/SRX322285_germline_snvs_indels.vcf \
--variant:SRX322286 normal_cd34_caller/SRX322286_germline_snvs_indels.vcf \
--variant:SRX322287 normal_cd34_caller/SRX322287_germline_snvs_indels.vcf \
--variant:SRX322288 normal_cd34_caller/SRX322288_germline_snvs_indels.vcf \
--variant:SRX322289 normal_cd34_caller/SRX322289_germline_snvs_indels.vcf \
--variant:SRX322290 normal_cd34_caller/SRX322290_germline_snvs_indels.vcf \
--variant:SRX322291 normal_cd34_caller/SRX322291_germline_snvs_indels.vcf \
--variant:SRX322292 normal_cd34_caller/SRX322292_germline_snvs_indels.vcf \
--variant:SRX322293 normal_cd34_caller/SRX322293_germline_snvs_indels.vcf \
--variant:SRX322294 normal_cd34_caller/SRX322294_germline_snvs_indels.vcf \
--variant:SRX322295 normal_cd34_caller/SRX322295_germline_snvs_indels.vcf \
--variant:SRX322296 normal_cd34_caller/SRX322296_germline_snvs_indels.vcf \
--variant:SRX322297 normal_cd34_caller/SRX322297_germline_snvs_indels.vcf \
--variant:SRX322298 normal_cd34_caller/SRX322298_germline_snvs_indels.vcf \
-o normal_cd34_caller/PON_caller_target_regions.vcf \
-L ./FindingSuitableCoverageData/Leucegene_target_genes.bed \ # restrict to regions of interest
-genotypeMergeOptions UNIQUIFY 

gatk  -T VariantsToTable \
-R path_hg19_reference/genome.fa \
-V normal_cd34_caller/PON_caller_target_regions.vcf \
--splitMultiAllelic \
-F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F set -GF AF \
-o normal_cd34_caller/PON_caller_target_regions.table
```

The code below was used to create the final PON for each pipeline. The `parsePON` function belong to the package `samplepower` which was developed to analyse the output from every downsampled run in this analysis and to compute the analysis for the paper. `samplepower` can be installed with the following code:

``` r
# https://github.com/annaquaglieri16/samplepower
library(devtools)
devtools::install_github("annaquaglieri16/samplepower")
library(samplepower)
```

The `caller` arguments is one of `varscan`, `mutect` or `vardict`

``` r
# VarScan2
var_linkPON <- file.path("PON_varscan_target_regions.table")
var_pon <- samplepower::parsePON(var_linkPON,caller="varscan")
var_pon$Location <- paste(var_pon$CHROM,var_pon$POS,sep="_")
```

### Annotation databases and genomic features

#### COSMIC, dbSNP and ExAC

The VEP 0.89 was used to annotate the variants output from the callers. VEP comes with a cache folders containing annotations associated with each genome build. We used VEP annotations to define if a variant was present in the COSMIC, dbSNP and ExAC databases.

#### RNA editing sites

We downloaded the database of RNA editing sites for the hg19 human reference genome from <http://rnaedit.com/download/> and imported and processed in `R` with the following code. The `variants` dataframe is saved in `FindingSuitableCoverageData/RADAR_hg19_allRNAedit.txt`.

``` r
library(GenomicRanges)

RADAR_hg19_allRNAedit <- read_table2("FindingSuitableCoverageData/RADAR_hg19_allRNAedit.txt", col_names = FALSE)
colnames(RADAR_hg19_allRNAedit) <- c("chrom","start","end")

RADAR_hg19_allRNAedit <- subset(RADAR_hg19_allRNAedit,chrom %in% chroms)
RADAR_hg19_allRNAedit <- subset(RADAR_hg19_allRNAedit,!is.na(RADAR_hg19_allRNAedit$start) & !is.na(RADAR_hg19_allRNAedit$end))
RNAeditGR <- GRanges(RADAR_hg19_allRNAedit)
```

#### RepeatMasker

We downloaded the RepeatMasker track from the UCSC Genome Browser for the hg19 human reference genome and saved as the file `FindingSuitableCoverageData/hg19_repeat_regions.bed`.

``` r
# Repetitive regions
hg19_repeat_regions <- read_table2("FindingSuitableCoverageData/hg19_repeat_regions.bed", col_names = FALSE)
colnames(hg19_repeat_regions) <- c("chrom","start","end","name","X5","strand")
hg19_repeat_regions <- subset(hg19_repeat_regions,chrom %in% chroms)
repeatsGR <- GRanges(hg19_repeat_regions)
```

#### Exon boundaries

The range object containing the 4bp flanking exon boundaries was created with the code below and saved in `FindingSuitableCoverageData/hg19_exons_combined.rds`. The NCBI annotation for the gene Symbols was downloaded from `ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/` and saved in `FindingSuitableCoverageData/Homo_sapiens.gene_info_29-11-2017.gz`.

``` r
library(GenomicRanges)
library(Rsubread)

ncbi <- read.delim(gzfile(file.path("FindingSuitableCoverageData/Homo_sapiens.gene_info_29-11-2017.gz")))

hg19 <- Rsubread::getInBuiltAnnotation(annotation="hg19")   
chroms <- c(paste0("chr",1:22),"chrX","chrY","chrM")
match_symbol <- match(hg19$GeneID,ncbi$GeneID)
hg19$Symbol <- ncbi$Symbol[match_symbol]
hg19_exons_start <- GRanges(seqnames = hg19$Chr,IRanges(start=hg19$Start - 4,end=hg19$Start),strand = hg19$Strand,
                            mcols=data.frame(side=rep("left",nrow(hg19))))
hg19_exons_end <- GRanges(seqnames = hg19$Chr,IRanges(start=hg19$End,end=hg19$End+4),strand = hg19$Strand,
                          mcols=data.frame(side=rep("right",nrow(hg19))))
hg19_exons_combined <- append(hg19_exons_start,hg19_exons_end)
```

#### Homopolymers

The code below was used to create the ranges of homopolymers. The hg19 UCSC reference genome FASTA file was downloaded from `http://sapac.support.illumina.com/sequencing/sequencing_software/igenome.html`.

``` r
library(seqinr)
chroms <- c(paste0("chr",1:22),"chrX","chrY","chrM")

fastaHg19 <- read.fasta("ref_genome.fa")

fastaHg19_chr <- fastaHg19[names(fastaHg19) %in% chroms]

# Make an RLe for every chromosome
rle_Hg19 <- lapply(fastaHg19_chr,function(chrom){
  Rle(as.numeric(as.factor(chrom)))
})

# for every element of an RLE for every chromosome save the start, the end, and the length of the rle run
rle_Hg19_homop <- lapply(rle_Hg19,function(chrom){
  end = end(chrom)
  start = start(chrom)
  len <- runLength(chrom)
  data.frame(start,end,len) %>%
    filter(len > 5)
})

# Combine them all to create a GRanges object
rle_chrom <- do.call(rbind,rle_Hg19_homop)
rle_chrom$chrom <- rep(names(rle_Hg19_homop),
                       times=sapply(rle_Hg19_homop,nrow))

GRanges_homop <- GRanges(seqnames=rle_chrom$chrom,
                         IRanges(start=rle_chrom$start,
                                 end=rle_chrom$end))
```

Bibliography
============

Lavallée, Vincent-Philippe, Sébastien Lemieux, Geneviève Boucher, Patrick Gendron, Isabel Boivin, Richard Neil Armstrong, Guy Sauvageau, and Josée Hébert. 2016. “RNA-sequencing Analysis of Core Binding Factor AML Identifies Recurrent Zbtb7a Mutations and Defines Runx1-Cbfa2t3 Fusion Signature.” *Blood*, March.
