#!/bin/bash

# This script uses hisat2 to align single-end reads.
# It sorts the mapped reads by name and coordinates.
# instructions on https://github.com/cyaguesa/SL-quant

# REQUIREMENTS (see installation instruction if unmet):      VERSION USED 
# - samtools installed and in your path                      1.5
# - trimmomatic installed and in your path                   0.36        
# - hisat2 installed and in your path                        2.0.5             


# by Carlo Yague-Sanz, 2017


# PARAMETERS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

index="data/ce10_hisat2_index/genome"  # genome index file.
core=1                                 # number of cores used for sorting.
mem="3328M"                            # amount of memory by core used for sorting.
set -e                                 # stops script at the first error.
threads=4                              # number of threads to use for mapping.

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

echo ""

if [ $# -eq 0 ]; then
  echo "   Usage : ./map_reads_modENCODE.sh [file1.fq.gz] [[file2.fq.gz, ...]] "; echo ""
  echo "   Example : ./map_reads_modENCODE.sh data/reads/modENCODE/modENCODE_4594.fastq.gz"
  exit 1

else
  for i in $*
    do
    base=$(echo ${i} | rev | cut -c 10- | rev)
    echo "  fastq.gz input file: $base"
    mkdir -p $base

    echo "  trim adapters"
    trimmomatic SE $base.fastq.gz ${base}/trimmed.fq.gz ILLUMINACLIP:data/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "  map SE reads"
    hisat2 -p $threads --no-softclip --min-intronlen 20 --max-intronlen 5000 --rna-strandness R -x $index -U ${base}/trimmed.fq.gz -S ${base}/mapped.sam
    samtools view -b -F 260 ${base}/mapped.sam> ${base}/accepted_hits.bam
    samtools view -b -f 4 ${base}/mapped.sam> ${base}/unmapped.bam
    
    echo "... done... sorting..."
    samtools sort -@ $core -m $mem ${base}/accepted_hits.bam > ${base}/accepted_hits_sorted.bam

    echo "... done... indexing and cleaning..."
    samtools index ${base}/accepted_hits_sorted.bam
    samtools flagstat ${base}/accepted_hits_sorted.bam > ${base}/accepted_hits.bam.flagstat
    rm ${base}/accepted_hits.bam
    rm ${base}/mapped.sam

  done
fi
echo ""
