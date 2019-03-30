#!/bin/bash

# ## Building a custom FUsionFilter Dataset
# # https://github.com/FusionFilter/FusionFilter/wiki/Building-a-Custom-FusionFilter-Dataset

# datasets
gtf=$1
genome_fasta=$2
star_fusion_data=$3


mkdir -p ${star_fusion_data}

cd ${star_fusion_data}

# 0. Download fusion annotation
# check if fusion annotation is present or get it from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
annotation=./Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v26_CTAT_lib_July192017.source_data/fusion_lib.dat.gz

if [[ -r ${annotation} ]]; then 

	echo "Fusion annotation has already been downloaded"

else


	wget -r -np --no-host-directories https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v26_CTAT_lib_July192017.source_data.tar.gz -P .

	tar xvzf ./Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v26_CTAT_lib_July192017.source_data.tar.gz 

	annotation=./Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v26_CTAT_lib_July192017.source_data/fusion_lib.dat.gz

fi

# 1. Extract cDNA sequences
/home/users/allstaff/quaglieri.a/software/STAR-Fusion/FusionFilter/util/gtf_file_to_cDNA_seqs.pl \
 ${gtf} ${genome_fasta} > ${star_fusion_data}/cDNA_seqs.fa 


# 3. All-vs-all BLASTN search
# make the cDNA_seqs.fa file blastable
chmod 775 *

makeblastdb -in ${star_fusion_data}/cDNA_seqs.fa -dbtype nucl  \
> ${star_fusion_data}/log_makeblastdb.txt

# perform the blastn search
chmod 775  *

## Without RepeatMasker output
blastn -query ${star_fusion_data}/cDNA_seqs.fa \
-db ${star_fusion_data}/cDNA_seqs.fa \
-max_target_seqs 10000 -outfmt 6 \
-evalue 1e-3 -lcase_masking \
-num_threads 10 \
-word_size 11  >  ${star_fusion_data}/blast_pairs.outfmt6

chmod 775 * 


#####

# Replace the transcript identifiers in the blast output with gene symbols (and gzipping output) by then running:
~/software/STAR-Fusion/FusionFilter/util/blast_outfmt6_replace_trans_id_w_gene_symbol.pl \
${star_fusion_data}/cDNA_seqs.fa \
${star_fusion_data}/blast_pairs.outfmt6  | \
gzip > ${star_fusion_data}/blast_pairs.gene_syms.outfmt6.gz 

# # Prep the Custom FusionFilter Dataset
# ## Need to update sometihng
# module add gmap-gsnap/2016-04-04

# Old version
# cd ${star_fusion_data}
# /home/users/allstaff/quaglieri.a/software/STAR-Fusion/FusionFilter/prep_genome_lib.pl \
# --genome_fa ${genome_fasta} \
# --gtf ${gtf} \
# --blast_pairs ./blast_pairs.gene_syms.outfmt6.gz \
# --cdna_fa ./cDNA_seqs.fa \
# --CPU 10 

# Create star-fusion with annotation 401
cd ${star_fusion_data}
/home/users/allstaff/quaglieri.a/software/STAR-Fusion/FusionFilter/prep_genome_lib.pl \
--genome_fa ${genome_fasta} \
--gtf ${gtf} \
--blast_pairs ./blast_pairs.gene_syms.outfmt6.gz \
--fusion_annot_lib ${annotation} \
--CPU 10 

/home/users/allstaff/quaglieri.a/software/STAR-Fusion/FusionFilter/util/index_pfam_domain_info.pl  \
--pfam_domains PFAM.domtblout.dat.gz \
--genome_lib_dir ctat_genome_lib_build_dir

# done
