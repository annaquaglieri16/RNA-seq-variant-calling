#!/bin/bash

module add STAR/2.5 

genome=$1
gtf=$2
genomeout=$3
genome_version=$4
maxSpan=$5

set -e

genomedir=${genomeout}/star_index_${genome_version}_${maxSpan}

mkdir -p ${genomedir}

STAR \
--runThreadN $(nproc) \
--runMode genomeGenerate \
--genomeDir ${genomedir} \
--genomeFastaFiles ${genome} \
--sjdbOverhang ${maxSpan} \
--sjdbGTFfile ${gtf}


