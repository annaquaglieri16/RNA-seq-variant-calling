## Check read groups in the original file

fastqdir=/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/quaglieri.a/fastq
find . -name "*-*.fastq" | cut -f4 -d "_" | sort | uniq | sed s/-/+/ > barcodes.txt

for combined in ${fastqdir}/*combined*.fastq; do

	out=$(basename $combined .fastq)
	grep -f barcodes.txt $combined | cut -f2 -d$' ' | cut -f4 -d ":" >> ${fastqdir}/${out}_read_names.txt

	read_names=${fastqdir}/${out}_read_names.txt
	echo "$read_names" >> ${fastqdir}/summary_readnames.txt
	cat $read_names | sort | uniq >> ${fastqdir}/summary_readnames.txt

done 



for fastq in $(find . -name "*.fastq" | sort);do

	nreads=$(wc -l $fastq)
	echo "$nreads"
	echo "$nreads" >>  ${fastqdir}/nreads.txt

done
# read numbers are all good!


## Possible problems - no problems
6219130_ENTRY_combined_R2
GAGATTCCGAGATTCC+GGCTCTGA
GAGATTCC+GGCTCTGA


## 
grep AGCGATAG+CCTATAGCGATAG+CCTATCCT 1152241_ENTRY_combined_R1.fastq
grep AGCGATAG+CCTAGCGATAG+CCTATCCT 1152241_ENTRY_combined_R1.fastq
grep GAGATTCCGAGATTCC+GGCTCTGA 6219130_ENTRY_combined_R2.fastq

# nothing 

