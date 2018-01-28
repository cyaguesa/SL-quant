#!/bin/bash

# This script uses tophat2 to align paired-end reads.
# It sorts the mapped reads by name and coordinates.
# instructions on https://github.com/cyaguesa/SL-quant

# REQUIREMENTS (see installation instruction if unmet):      VERSION USED       NOTE
# - TopHat2 installed and in your path                       2.1.1
# - samtools installed and in your path                      1.5
# - trimmomatic installed with path set in parameters        0.36               Manually set the path.


# by Carlo Yague-Sanz, 2017


# PARAMETERS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

index="data/ce10_bowtie2_index/genome" # genome index file
core=1                                 # number of cores used for sorting
mem="3328M"                            # amount of memory by core used for sorting
set -e                                 # stops script at the first error
trimmomatic="/Users/Carlo/Documents/binaries/Trimmomatic-0.36/trimmomatic-0.36.jar"           # path to trimmomatic                           

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
    java -jar $trimmomatic SE $base.fastq.gz ${base}/trimmed.fq.gz ILLUMINACLIP:data/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "  map SE reads"
    tophat2 -r 200 --no-sort -N 2 --read-edit-dist 2 --min-intron-length 5 --max-intron-length 3000 --max-multihits 4 --no-discordant --no-coverage-search -p 4 --library-type fr-firststrand --output-dir $base $index ${base}/trimmed.fq.gz

    echo "... done... sorting..."
    samtools sort -@ $core -m $mem ${base}/accepted_hits.bam > ${base}/accepted_hits_sorted.bam

    echo "... done... indexing and cleaning..."
    samtools index ${base}/accepted_hits_sorted.bam
    samtools flagstat ${base}/accepted_hits_sorted.bam > ${base}/accepted_hits.bam.flagstat
    rm ${base}/accepted_hits.bam

  done
fi
echo ""