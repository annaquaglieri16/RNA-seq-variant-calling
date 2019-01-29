## freebayes

freebayes --bam bamfile --vcf -f genome_fasta --targets regions_bed --min-base-quality 20 --min-mapping-quality 20 \
--min-alternate-count 2 --min-coverage 10


# varscan
--min-coverage 10 --min-reads2 2 --variants 1 --output-vcf 1 --min-avg-qual 20 --min-var-freq 0.01

# vardict
AF_THR="0.01" # minimum allele frequency
O - average min mapping quality for reads
# regions should be provided as chr:start-end
# -r 2 minimum n variant reads
VarDict -G genome_fasta -f 0.01 -N sample_name -b bamfile -c 1 -S 2 -E 3 -g 4 -O 20 -q 20 -R regions -r 2 -t -th 10 -v | teststrandbias.R | var2vcf_valid.pl -N sample_name -E -f 0.01

AF_THR="0.01" # minimum allele frequency
vardict -G genome_fasta -f 0.01 -N tumor_sample_name -b "tumour.bam|normal.bam" -c 1 -S 2 -E 3 -g 4 -O 20 -q 20 -R regions -r 2 -t -th 10 -v | testsomatic.R | var2vcf_paired.pl -N "tumor_sample_name|normal_sample_name" -f $0.01

