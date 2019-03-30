# cat tmp_header tmp_bottom_ALT > SRX322282_mpileup2cns_snvs_indels_ALT.vcf

gatk -T CombineVariants \
-R $genome_fasta \
--variant:SRX322282 SRX322282_mpileup2cns_snvs_indels.vcf \
--variant:SRX322283 SRX322283_mpileup2cns_snvs_indels.vcf \
--variant:SRX322284 SRX322284_mpileup2cns_snvs_indels.vcf \
--variant:SRX322285 SRX322285_mpileup2cns_snvs_indels.vcf \
--variant:SRX322286 SRX322286_mpileup2cns_snvs_indels.vcf \
--variant:SRX322287 SRX322287_mpileup2cns_snvs_indels.vcf \
--variant:SRX322288 SRX322288_mpileup2cns_snvs_indels.vcf \
--variant:SRX322289 SRX322289_mpileup2cns_snvs_indels.vcf \
--variant:SRX322290 SRX322290_mpileup2cns_snvs_indels.vcf \
--variant:SRX322291 SRX322291_mpileup2cns_snvs_indels.vcf \
--variant:SRX322292 SRX322292_mpileup2cns_snvs_indels.vcf \
--variant:SRX322293 SRX322293_mpileup2cns_snvs_indels.vcf \
--variant:SRX322294 SRX322294_mpileup2cns_snvs_indels.vcf \
--variant:SRX322295 SRX322295_mpileup2cns_snvs_indels.vcf \
--variant:SRX322296 SRX322296_mpileup2cns_snvs_indels.vcf \
--variant:SRX322297 SRX322297_mpileup2cns_snvs_indels.vcf \
--variant:SRX322298 SRX322298_mpileup2cns_snvs_indels.vcf \
-o combined_variants_control_cd341.vcf \
-genotypeMergeOptions UNIQUIFY 


# Filter normals when they have no variant allele
normal_bam=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/hg19_aligned/variant_calling/samtools/combined_variants_control_cd341.vcf
awk '($0 !~ /^#/)' $normal_bam > without_header_normals.vcf
awk '{ if ($5 != ".") print $0 }' without_header_normals.vcf > filtered_normals_no_header.vcf
head -55 $normal_bam > header_normals.vcf
cat header_normals.vcf filtered_normals_no_header.vcf > combined_variants_control_cd341_filtered.vcf


# process variants
variantdir=${star_2pass}/../variant_calling/samtools

${workdir}/process_variants_all_genes.sh ${star_2pass}/../variant_calling/samtools/combined_variants_control_cd341_filtered.vcf "PON_CD34"


# process combined variants
vcf=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/hg19_aligned/variant_calling/samtools/combined_variants_control_cd341.vcf
${workdir}/process_variants.sh $vcf "combined_variants_control_cd341"

awk '($5 != ".")' combined_variants_control_cd341_variants_annotated.table > combined_variants_control_cd341_variants_annotated_ALT.table

awk '(print $1-$15)' combined_variants_control_cd341_variants_annotated_ALT.table

cut combined_variants_control_cd341_variants_annotated_ALT.table -f 1-15 > \
combined_variants_control_cd341_variants_annotated_ALT_noSamples.table

###
tables=$(find . -maxdepth 1 -name "*_variants_annotated.table" | sort)
for table in $tables; do

out=${table/_variants_annotated.table/_variants_annotated_ALT.table}

awk '($5 != ".")' $table > $out

done

echo "Control initial done" | mail -s "Control initial done" quaglieri.a@wehi.edu.au

