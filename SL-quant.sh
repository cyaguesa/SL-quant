#!/bin/bash

# In single-end mode, SL-quant identifies trans-splicing events as unmapped reads containing a splice leader (SL) sequence at the 5’ end of the read. 
# After trimming of those SL sequences, the reads are re-mapped on the genome and counted at the gene level.
# In paired-end mode, SL-quant identifies trans-splicing events as read pairs with one end unmapped and starting with a SL sequence. 
# Using the positional information of the mapped ends of the pairs, it outputs the number of trans-splicing events per genes.
# instructions on https://github.com/cyaguesa/SL-quant

# REQUIREMENTS (see installation instruction if unmet):      VERSION USED       NOTE

# - blastn (blast+) installed and in your path               2.6.0+
# - samtools installed and in your path                      1.5
# - picard installed and in your path                        2.9.0
# - bowtie2 installed and in your path                       2.3.2              Only required on single-end mode.
# - featureCounts installed and in your path                 1.5.0              


# by Carlo Yague-Sanz, 2018


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# PARAMETERS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

set -e                                 # stops script at the first error.
SINGLE="single"                        # set to "single" for single-end mode. Any other value for paired-end mode.
SL_db="data/blast_db/SL.fasta"         # path to SL sequence database (for blast).
gene_annotation="data/genes.SAF"       # annotation file
index="data/ce10_bowtie2_index/genome" # genome index file (only required on single-end mode).
paired_orientation="fr-firststrand"    # ignored in single-end mode. Value={"fr-firststrand" (default), "fr-secondstrand", "fr-unstrand"}
single_orientation="stranded"          # ignored in paired-end mode. Value={"stranded" (default), "reversely_stranded", "unstranded"}


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# PAIRED-END MODE
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

if [ "$SINGLE" != "single" ]; then
  if [ $# != 3 ]; then
     echo ""; echo "   Usage : ./SL-quant.sh [mapped.bam] [unmapped.bam] [ouput_dir/base]"; echo ""
     echo "   Example : ./SL-quant.sh data/reads/SRR1585277/accepted_hits_sorted.bam data/reads/SRR1585277/unmapped.bam SL-quant_results/SRR1585277";echo""
     exit 1

  else
    echo ""
    echo "SL-quant.sh: PAIRED-END MODE"
    echo "   [mapped.bam] = $1"
    echo "   [unmapped.bam] = $2"
    echo "   [ouput_dir/base] = $3"
    echo "   [SL blast database] = $SL_db"
    echo "   [gene annotation] = $gene_annotation"
    echo "   [read orientation] = $paired_orientation";echo ""

    echo "   fetch read pairs with one end unmapped..."
    
    if [ "$paired_orientation" == "fr-firststrand" ]; then
 
      echo "      [1/2] get R2 reads unmapped with mate mapped and convert to fasta"
      samtools view -f 133 -F 8 ${2} | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${3}_oneEnd_unmapped.fasta

      echo "      [2/2] get primary alignments of R1 reads mapped with mate unmapped"
      samtools view -u -f 73 -F 260 ${1}  > ${3}_oneEndMapped.bam

      featureCounts_S=2   # set R1 read orientation for featureCounts (reverse)

    elif [ "$paired_orientation" == "fr-secondstrand" ]; then

      echo "      [1/2] get R1 reads unmapped with mate mapped and convert to fasta"
      samtools view -f 69 -F 8 ${2} | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${3}_oneEnd_unmapped.fasta

      echo "      [2/2] get primary alignments of R2 reads mapped with mate unmapped"
      samtools view -u -f 137 -F 260 ${1}  > ${3}_oneEndMapped.bam

      featureCounts_S=2   # set R2 read orientation for featureCounts (reverse)

    else

      echo "      [1/2] get reads unmapped with mate mapped and convert to fasta"
      samtools view -f 5 -F 8 ${2} | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${3}_oneEnd_unmapped.fasta

      echo "      [2/2] get reads mapped with mate unmapped"
      samtools view -u -f 9 -F 260 ${1}  > ${3}_oneEndMapped.bam

      featureCounts_S=0   # set read orientation for featureCounts (unstranded)

    fi

    echo "   done... blast unmapped reads on SL sequences..."
    blastn -query ${3}_oneEnd_unmapped.fasta -task blastn -db $SL_db  -outfmt 6 -max_target_seqs 1 -num_threads 4 -word_size 8 > ${3}_blasted.txt 2>${3}_log.txt

    echo "   done... filter significant alignment with qstart==1..."
    awk '$7 == 1 && $11 < "0.05" {print $0}' ${3}_blasted.txt | grep "SL1" > ${3}_blasted_SL1.txt
    awk '$7 == 1 && $11 < "0.05" {print $0}' ${3}_blasted.txt | grep "SL2" > ${3}_blasted_SL2.txt

    echo "   done... retrieve mapped mates..."

    echo "      [1/2] of SL1-containing reads"
    cut -f 1 ${3}_blasted_SL1.txt > ${3}_SL1_IDs.txt
    picard FilterSamReads FILTER=includeReadList READ_LIST_FILE=${3}_SL1_IDs.txt I=${3}_oneEndMapped.bam O=${3}_SL1_mates.bam 2>> ${3}_log.txt

    echo "      [2/2] of SL2-containing reads"
    cut -f 1 ${3}_blasted_SL2.txt > ${3}_SL2_IDs.txt
    picard FilterSamReads FILTER=includeReadList READ_LIST_FILE=${3}_SL2_IDs.txt I=${3}_oneEndMapped.bam O=${3}_SL2_mates.bam 2>> ${3}_log.txt

    echo "   done... summarize..."
    featureCounts -s $featureCounts_S -F SAF -g GeneID -T 4 -a $gene_annotation -o ${3}_counts.txt ${3}_SL1_mates.bam ${3}_SL2_mates.bam 2>> ${3}_log.txt

  fi


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# SINGLE-END MODE
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

else
  if [ $# != 2 ]; then
    echo ""; echo "   Usage : ./SL-quant.sh [unmapped.bam] [ouput_dir/base]"; echo ""
    echo "   Example : ./SL-quant.sh data/reads/SRR1585277/unmapped.bam SL-quant_results/SRR1585277_single"
    echo "   Example2 : ./SL-quant.sh data/reads/modENCODE/modENCODE_4594/unmapped.bam SL-quant_results/modENCODE_4594"
    exit 1

  else
    echo ""
    echo "SL-quant.sh: SINGLE-END MODE"
    echo "   [unmapped.bam] = $1"
    echo "   [ouput_dir/base] = $2"
    echo "   [SL blast database] = $SL_db"
    echo "   [gene annotation] = $gene_annotation"
    echo "   [read orientation] = $single_orientation";echo ""

    paired_entries=$(samtools view -c -f 1 $1) # set this to 0 to save time if you don't want to check if the data is paired
    if [ $paired_entries -gt 0 ]; then
      echo ""; echo "WARNING : there are ${paired_entries} 'paired in sequencing' reads in ${1}. Consider running the script in paired-end mode."; echo ""
    fi

    echo "   blast unmapped reads on SL sequences..."
    samtools view -f 4 ${1} | awk '{OFS="\t"; print ">"$1"\n"$10}' | blastn -db $SL_db -outfmt 6 -max_target_seqs 1 -num_threads 4 -word_size 8 > ${2}_blasted.txt 2>${2}_log.txt

    echo "   done... filter..."
    awk '$7 == 1 && $11 < "0.05" {print $0}' ${2}_blasted.txt | grep "SL1" > ${2}_blasted_SL1.txt
    awk '$7 == 1 && $11 < "0.05" {print $0}' ${2}_blasted.txt | grep "SL2" > ${2}_blasted_SL2.txt

    echo "   done... retrieve SL-containing reads"

    echo "      [1/2] of SL1-containing reads"
    cut -f 1 ${2}_blasted_SL1.txt > ${2}_SL1_IDs.txt
    picard FilterSamReads FILTER=includeReadList VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 QUIET=TRUE READ_LIST_FILE=${2}_SL1_IDs.txt I=${1} O=${2}_SL1.sam


    echo "      [2/2] of SL2-containing reads"
    cut -f 1 ${2}_blasted_SL2.txt > ${2}_SL2_IDs.txt
    picard FilterSamReads FILTER=includeReadList VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 QUIET=TRUE READ_LIST_FILE=${2}_SL2_IDs.txt I=${1} O=${2}_SL2.sam

    echo "   done... trim SL-containing reads and convert to fastq"
    paste <(samtools view ${2}_SL1.sam | cut -f 1-11) <(cat ${2}_blasted_SL1.txt) | awk '{OFS="\t"; print $1,"4",$3,$4,$5,$6,$7,$8,$9, substr($10, $19),  substr($11, $19)}' | picard SamToFastq VALIDATION_STRINGENCY=SILENT QUIET=TRUE I=/dev/stdin FASTQ=${2}_SL1_trimmed.fq
    paste <(samtools view ${2}_SL2.sam | cut -f 1-11) <(cat ${2}_blasted_SL2.txt) | awk '{OFS="\t"; print $1,"4",$3,$4,$5,$6,$7,$8,$9, substr($10, $19),  substr($11, $19)}' | picard SamToFastq VALIDATION_STRINGENCY=SILENT QUIET=TRUE I=/dev/stdin FASTQ=${2}_SL2_trimmed.fq

    echo "   done... map with bowtie2"
    bowtie2 --no-unal --quiet -U ${2}_SL1_trimmed.fq -x $index | samtools view -b -F 256 > ${2}_SL1_remapped.bam
    bowtie2 --no-unal --quiet -U ${2}_SL2_trimmed.fq -x $index | samtools view -b -F 256 > ${2}_SL2_remapped.bam

    echo "   done... summarize..."
    if [ "$single_orientation" == "stranded" ]; then               # set read orientation for featureCounts (forward)
      featureCounts_S=1   
    elif [ "$single_orientation" == "reversely_stranded" ]; then   # set read orientation for featureCounts (reverse)
      featureCounts_S=2   
    else                                                           # set read orientation for featureCounts (unstranded)
      featureCounts_S=0
    fi
    featureCounts -s $featureCounts_S -F SAF -g GeneID -T 4 -a $gene_annotation -o ${2}_counts.txt ${2}_SL1_remapped.bam ${2}_SL2_remapped.bam 2>> ${2}_log.txt
  fi

fi
echo "   All done ! Have a nice day :)"
echo ""