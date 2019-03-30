#!/bin/bash

genome=$1
gtf=$2
genomeout=$3
maxSpan=$4

set -e

mkdir -p ${genomeout}

STAR \
--runThreadN $(nproc) \
--runMode genomeGenerate \
--genomeDir ${genomeout} \
--genomeFastaFiles ${genome} \
--sjdbOverhang ${maxSpan} \
--sjdbGTFfile ${gtf}


