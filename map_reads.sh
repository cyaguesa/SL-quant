#!/bin/bash

# This script uses hisat2 to align paired-end reads.
# It sorts the mapped reads by name and coordinates.
# instructions on https://github.com/cyaguesa/SL-quant

# REQUIREMENTS (see installation instruction if unmet):      VERSION USED 
# - hisat2 installed and in your path                        2.0.5             
# - samtools installed and in your path                      1.5

# by Carlo Yague-Sanz, 2017


# PARAMETERS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

dir="data/reads/SRR1585277"            # data directory
index="data/ce10_hisat2_index/genome"  # genome index file
offset=12                              # length of "_1.fastq.gz" +1. When the ".R1.fastq.gz" notation is used, this should be set to 13.
core=1                                 # number of cores used for sorting
mem="3328M"                            # amount of memory by core used for sorting
set -e                                 # stops script at the first error
threads=4                              # number of threads to use for mapping.

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––


for i in $(ls $dir/*.fastq.gz | rev | cut -c $offset- | rev | uniq)
  do
  echo ""; echo "fastq file : ${i}"

  echo "... aligning..."
  mkdir -p ${i}
  hisat2 -p $threads --no-softclip --no-discordant --min-intronlen 20 --max-intronlen 5000 --rna-strandness FR -x $index -1 ${i}_1.fastq.gz -2 ${i}_2.fastq.gz -S ${i}/mapped.sam
  samtools view -b -F 4 ${i}/mapped.sam> ${i}/accepted_hits.bam
  samtools view -b -f 4 ${i}/mapped.sam> ${i}/unmapped.bam

  echo "... done... sorting..."
  samtools sort -@ $core -m $mem ${i}/accepted_hits.bam > ${i}/accepted_hits_sorted.bam
#  samtools sort -n -@ $core -m $mem ${i}/accepted_hits.bam > ${i}/accepted_hits_Nsorted.bam

  echo "... done... indexing and cleaning..."
  samtools index ${i}/accepted_hits_sorted.bam
  samtools flagstat ${i}/accepted_hits_sorted.bam > ${i}/accepted_hits.bam.flagstat
  rm ${i}/accepted_hits.bam
  rm ${i}/mapped.sam

  echo "... done. Have a nice day !"
  done
