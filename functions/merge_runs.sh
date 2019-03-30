#!/bin/bash

################
# Merge bamfiles
################

module add sambamba

sample=$1 
bamdir=$2

set -e

merged_runs=${bamdir}_merged_runs
mkdir -p ${merged_runs}

merged_bam=${merged_runs}/${sample}merged_Dupl_rg.bam
merged_DeDupl_bam=${merged_runs}/${sample}merged_DeDupl_rg.bam

############################
# Bam with Marked Duplicates
############################

if [[ -s ${merged_bam} ]] ; then 

	echo "Merge Marked Dupl exists for $(basename ${merged_bam})"

else

	echo "merging runs for sample $sample"
	sambamba merge ${merged_runs}/${sample}merged_Dupl_rg.bam ${bamdir}/*${sample}*Aligned.reorderedDupl.rg.bam -t $(nproc)

	## Convert BAM to CRAM
	#sambamba view -T $genome_fasta -f cram -o ${star_2pass}/${sample}_merged_Dupl_rg.bam.cram ${star_2pass}/${sample}_merged_Dupl_rg.bam --nthreads $(nproc)
	samtools index ${merged_runs}/${sample}merged_Dupl_rg.bam

fi 

##################
# DeDuplicated bam
##################

if [[ -s ${merged_DeDupl_bam} ]] ; then 
 
  	echo "Merge DeDupl exists for  (basename ${merged_DeDupl_bam})"
 
  else
 
  	echo "merging deduplicated runs for sample $sample"
 
  	sambamba merge ${merged_runs}/${sample}merged_DeDupl_rg.bam ${bamdir}/${sample}*Aligned.reordered.DeDupl.rg.bam -t $(nproc)
 
  	samtools index ${merged_runs}/${sample}merged_DeDupl_rg.bam
 
fi 


