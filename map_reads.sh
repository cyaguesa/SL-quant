#!/bin/bash

# This script uses tophat2 to align paired-end reads.
# It sorts the mapped reads by name and coordinates.
# instructions on https://github.com/cyaguesa/SL-quant

# REQUIREMENTS (see installation instruction if unmet):      VERSION USED 
# - tophat2 installed and in your path                       2.1.1
# - samtools installed and in your path                      1.5

# by Carlo Yague-Sanz, 2017


# PARAMETERS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

dir="data/reads/SRR1585277"            # data directory
index="data/ce10_bowtie2_index/genome" # genome index file
offset=12                              # length of "_1.fastq.gz" +1. When the ".R1.fastq.gz" notation is used, this should be set to 13.
core=1                                 # number of cores used for sorting
mem="3328M"                            # amount of memory by core used for sorting
set -e                                 # stops script at the first error

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––


for i in $(ls $dir/*.fastq.gz | rev | cut -c $offset- | rev | uniq)
  do
  echo ""; echo "fastq file : ${i}"

  echo "... aligning..."
  tophat --rg-sample ${i} --rg-id ${i} -r 200 --no-sort --mate-std-dev 150 -N 2 --read-edit-dist 2 --min-intron-length 5 --max-intron-length 3000 --max-multihits 4 --no-discordant --no-coverage-search -p 4 --library-type fr-firststrand --output-dir ${i} $index ${i}_1.fastq.gz ${i}_2.fastq.gz

  echo "... done... sorting..."
  samtools sort -@ $core -m $mem ${i}/accepted_hits.bam > ${i}/accepted_hits_sorted.bam
#  samtools sort -n -@ $core -m $mem ${i}/accepted_hits.bam > ${i}/accepted_hits_Nsorted.bam

  echo "... done... indexing and cleaning..."
  samtools index ${i}/accepted_hits_sorted.bam
  samtools flagstat ${i}/accepted_hits_sorted.bam > ${i}/accepted_hits.bam.flagstat
  rm ${i}/accepted_hits.bam

  echo "... done. Have a nice day !"
  done