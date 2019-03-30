## Merge fastqfiles
fastqdir=$1
scripts=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/scripts/RNA-seq_variant_calling
#bamdir=/home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/hg19_aligned/STAR_1pass_preliminary
samples_runs=$2
srx=$3

cd $fastqdir

nice -12 python ${scripts}/combined_fastqs.py $samples_runs $srx

# Merge bamfiles to be able to set separate read groups
# set -e

# nice -9 cat $srx | while read srx_sample; do

# run=$(grep $srx_sample $samples_runs | cut -f2 -d ' ' | xargs -I {} bash -c 'find /home/users/allstaff/quaglieri.a/PHD_project/GEO_Leucegene_data/hg19_aligned/STAR_1pass_preliminary/ -name {}Aligned.sortedByCoord.rg.out.bam')

# echo "$run" > ${bamdir}/${srx_sample}_bamlist.txt

# samtools merge -f ${bamdir}/${srx_sample}_combined.sorted.rg.out.bam -b ${bamdir}/${srx_sample}_bamlist.txt

# done

