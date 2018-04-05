#!/bin/bash

# This script provides statistics about SL-sites (provided as bam files).
# It also checks whether the consensus sequence is "AG" or not
# and output the sequence around SL-sites in fasta.
# instructions on https://github.com/cyaguesa/SL-quant

# REQUIREMENTS (see installation instruction if unmet):      VERSION USED
# - bedtools installed and in your path                            2.26.0


# by Carlo Yague-Sanz, 2018


# PARAMETERS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

genome_info="data/chrom_summary"                      # genome index file
set -e                                                # stops script at the first error
genome_fasta="data/ce10_hisat2_index/genome.fa"       # genome index file

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

echo ""

if [ $# -eq 0 ]; then
  echo "   Usage : ./SL_sites.sh [file1.bam] [[file2.bam, ...]] "; echo ""
  echo "   Example : ./SL_sites.sh SL-quant_results/*/*_remapped.bam"
  exit 1

else
  for i in $*
    do
    echo "  find number of SL-transplicing events by sites (not by genes): ${i}"
    paired_entries=$(samtools view -c -f 1 ${i})
    if [ $paired_entries -gt 0 ]; then
      echo "  [paired-end data]"
      samtools view -b -f 128 ${i} > ${i}_R2.bam
      input="${i}_R2.bam"
    else
      echo "  [single-end data]"
      input="${i}"
    fi


    events=$(samtools flagstat ${input} | head -n 1 | cut -f 1 -d " ")
    echo "        number of SL-transplicing events detected: ${events}"

    samtools sort ${input} | bedtools genomecov -strand + -5 -dz -ibam stdin -g ${genome_info} > ${i}_fwd.tab
    samtools sort ${input} | bedtools genomecov -strand - -5 -dz -ibam stdin -g ${genome_info} > ${i}_rev.tab
    awk 'BEGIN {OFS="\t"};{ print $1, $2-6, $2+5, ".", $3, "+" }' ${i}_fwd.tab > ${i}_fwd.bed
    awk 'BEGIN {OFS="\t"};{ print $1, $2-4, $2+7, ".", $3, "-" }' ${i}_rev.tab > ${i}_rev.bed
    cat ${i}_fwd.bed ${i}_rev.bed > ${i}.bed
    sites=$(wc -l < ${i}.bed)
    echo "        number of sites detected: ${sites}"

    bedtools getfasta -s -fi ${genome_fasta} -bed ${i}.bed > ${i}.fasta
    AG=$(awk '/^>/ {N++;next;} {print($0);}' ${i}.fasta | cut -c 5,6 | sort | uniq -c | grep -i "AG" | awk '{sum += $1} END {print sum}')
    ratio=$(echo "scale=3 ; 100 * $AG / $sites" | bc)
    echo "        number of sites which are 'AG' splice-sites: ${AG} (${ratio} %)"
    echo ""
    done
fi
echo ""

