#!/bin/bash

bamin=$1
sample=$2
bamdir=$(dirname $bamin)
sampleID=$3
# ID which marks the sample. Needed for calling variants with Mutect2 on bamfiles constituted of
# different read groups
set -u

mkdir -p ${bamdir}/ValidateSam/

# Check if validationSam has already been run and if the last file produced by this code exists
validSam=${bamdir}/ValidateSam/${sample}_validate.txt
bamRG=${bamdir}/${sample}_Aligned.sortedByCoord.rg.bam
bamDupl=${bamdir}/${sample}_Aligned.reorderedDupl.rg.bam
bamDuplIndex=${bamDupl}.bai
validSamDeDupl=${bamdir}/ValidateSam/${sample}_deDupl_validate.txt
bamDeDupl=${bamdir}/${sample}_Aligned.reordered.DeDupl.rg.bam
bamDeDuplIndex=${bamDeDupl}.bai

# If the validation has already been performed then the function won't remark duplicates or aggRG or validate the file

if [ -s $validSam ] ; then

	echo "Validation has already been performed for $sample see $validSam"

	# check the Validation is good
	cat $validSam | while read line; do

    	if [ "$line" == "No errors found" ] ; then

            echo "$sample is a valid bamfile"

            # Check if the file with marked duplicates exists
                
            if [ -s $bamDupl -a -s $bamDuplIndex ] ; then 

            	echo "$(basename $bamDupl) and $(basename $bamDuplIndex) exist and are bigger than 0."

            	else

            	echo "$(basename $bamDupl) or $(basename $bamDuplIndex) have been removed"

            fi

        else

            echo "Warning: the aligned $sample might be corrupted"

            exit
        
        fi

    done

else

	now=$(date)
	echo "$now" >> ${bamdir}/${sample}_pre_process_log


	if [ -s ${bamRG} ] ; then 

		echo "${bamRG} already exixts and it is > 0"

	else


		################################
		echo "Adding Reag Group to $sample" 
		echo "Adding Reag Group to $sample" >> ${bamdir}/${sample}_pre_process_log

		AddOrReplaceReadGroups \
		I=${bamin} \
		O=${bamRG} \
		RGID=${sample} \
		RGPU=${sample} \
		RGLB=${sampleID} \
		RGPL="illumina" \
		RGSM=${sampleID} 2> ${bamdir}/${sample}_addRG_Log

		sambamba index $bamRG

	fi

	if [ -s ${bamDupl} ] ; then 

		echo "${bamDupl} already exixts and it is > 0"

	else


		################################
	    echo "Marking duplicates to $sample"

		sambamba markdup ${bamRG} ${bamDupl} -t $(nproc) --show-progress 

		sambamba index ${bamDupl}

		echo "Duplicates flagged and index created for $sample"

	fi


	if [ -s ${bamDeDupl} ] ; then 

		echo "${bamDeDupl} already exixts and it is > 0"

	else


		################################
	    echo "Removing duplicates to $sample"

		sambamba markdup ${bamRG} ${bamDeDupl} -t $(nproc) --remove-duplicates --show-progress 

		sambamba index ${bamDeDupl}

		echo "Duplicates removed and index created for $sample"

	fi

	################################

	ValidateSamFile \
	INPUT=${bamDeDupl} \
	MODE=SUMMARY \
	IGNORE_WARNINGS=true \
	VALIDATE_INDEX=true \
	OUTPUT=${bamdir}/ValidateSam/${sample}_deDupl_validate.txt 

	ValidateSamFile \
	INPUT=${bamDupl} \
	MODE=SUMMARY \
	IGNORE_WARNINGS=true \
	VALIDATE_INDEX=true \
	OUTPUT=${bamdir}/ValidateSam/${sample}_validate.txt 

	echo "Validated $sample"
	echo "Validated $sample" >> ${bamdir}/${sample}_pre_process_log

fi



