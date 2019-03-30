

control_dir <- ""

 java -jar $GATK \
 -T CombineVariants \
 --arg_file inputs.list \
 -minN 2 \
 --setKey "null" \
 --filteredAreUncalled \
 --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
 -o 2_pon_combinevariants.vcf.gz \
 -R ~/Documents/ref/hg38/Homo_sapiens_assembly38.fasta \
 --genotypemergeoption UNIQUIFY