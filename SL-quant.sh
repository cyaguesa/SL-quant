#!/bin/bash

# In single-end mode, SL-quant identifies trans-splicing events as unmapped reads containing a splice leader (SL) sequence at the 5’ end of the read. 
# In paired-end mode, SL-quant identifies trans-splicing events as read pairs with one end unmapped and starting with a SL sequence while the other end is mapped. 
# After trimming of those SL sequences, the reads/pairs are re-mapped on the genome and counted at the gene level.
# instructions on https://github.com/cyaguesa/SL-quant

# REQUIREMENTS:	            VERSION USED 

# - blastn                    2.6.0
# - samtools                    1.5
# - picard                    2.9.0
# - featureCounts             1.5.0              
# - hisat2                    2.0.5             
# - cutadapt                   1.14              
# - bedtools                 2.26.0              


# by Carlo Yague-Sanz, 2017

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# USAGE
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

usage=$(cat <<EOF

SL-quant identifies trans-splicing events as unmapped reads containing a splice leader (SL) sequence at the 5’ end of the read. In paired-end mode, the unmapped reads are pre-filtered based on the status of their mate (it must be mapped). In sensitive mode, the criteria to identify SL sequences are less stringent, increasing sensitivity.

Detailed instructions on github.com/cyaguesa/SL-quant

USAGE: ./SL-quant.sh [-p -m mapped.bam] [-s] unmapped.bam output_base

Required arguments:
  unmapped.bam
				file containing unmapped reads.
  output_base
				base name (+path) for the ouput.

Optional arguments:
  --mapped mapped.bam, -m mapped.bam
				file containing mapped reads.
  --paired, -p
				run SL-quant in paired-end mode. 
				Requires -m argument.
  --sensitive, -s
				run SL-quant in sensitive mode.
  --help, -h
				show this help message and exit.


EOF)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# PARAMETERS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––



set -e                                 # stops script at the first error.
SL_db="data/blast_db/SL.fasta"         # path to SL sequence database (for blast).
gene_annotation="data/genes.gtf"       # annotation file
index="data/ce10_hisat2_index/genome"  # genome index file (only required on single-end mode).
paired_orientation="FR"                # ignored in single-end mode. Value={"FR" (default), "RF", "unstranded"}
single_orientation="R"                 # ignored in paired-end mode. Value={"F" (stranded), "R" (reversely stranded), "unstranded"}
threads=4                              # number of threads to use.
send=22                                # End of alignment in 'subject' threshold for blast.
align_length=5                         # length of the cutadapt alignment (default = 5). (paired-end mode only)


PARAMS=""
SINGLE="single"                        # set to "single" for single-end mode. Any other value for paired-end mode.
METHOD="specific"                      # method used to detect SL-containing reads. "specific" (default, with blast) or "sensitive" (with cutadapt)

while (( "$#" )); do
  case "$1" in
    -p|--paired)
      SINGLE="paired"
      shift
      ;;
    -m|--mapped)
      PARAMS="$PARAMS $2"
      shift 2
      ;;
    -s|--sensitive)
      METHOD="sensitive"
      shift
      ;;
    -h|--help)
      echo "$usage"
      exit 1
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      echo "$usage"
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

eval set -- "$PARAMS"


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# PAIRED-END MODE
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

if [ "$SINGLE" != "single" ]; then
  if [ $# != 3 ]; then
     echo "${usage}"
     exit 1

  else
    echo ""
    echo "SL-quant.sh: PAIRED-END MODE [$METHOD]"
    echo "   [mapped.bam] = $1"
    echo "   [unmapped.bam] = $2"
    echo "   [ouput_dir/base] = $3"
    echo "   [SL blast database] = $SL_db"
    echo "   [gene annotation] = $gene_annotation"
    echo "   [read orientation] = $paired_orientation";echo ""

    echo "   fetch read pairs with one end unmapped..."
    
    if [ "$paired_orientation" == "FR" ]; then
 
      echo "      [1/2] get unmapped R2 reads with mate mapped and convert to fastq"
      samtools view -u -f 133 -F 8 ${2} > ${3}_oneEnd_unmapped.bam

      echo "      [2/2] get primary alignments of mapped R1 reads with mate unmapped"
      samtools view -u -f 73 -F 260 ${1}  > ${3}_oneEndMapped.bam

      featureCounts_S=2                    # set R1 read orientation for featureCounts (reverse)
      strand="plus"                        # set blast database strandness (normal)

    elif [ "$paired_orientation" == "RF" ]; then

      echo "      [1/2] get R1 reads unmapped with mate mapped and convert to fastq"
      samtools view -f 69 -F 8 ${2} | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${3}_oneEnd_unmapped.fasta

      echo "      [2/2] get primary alignments of R2 reads mapped with mate unmapped"
      samtools view -u -f 137 -F 260 ${1}  > ${3}_oneEndMapped.bam

      featureCounts_S=2                    # set R2 read orientation for featureCounts (reverse)
      strand="minus"                       # set blast database strandness (reverse)

    else

      echo "      [1/2] get reads unmapped with mate mapped and convert to fastq"
      samtools view -f 5 -F 8 ${2} | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${3}_oneEnd_unmapped.fasta

      echo "      [2/2] get reads mapped with mate unmapped"
      samtools view -u -f 9 -F 260 ${1}  > ${3}_oneEndMapped.bam

      featureCounts_S=0                   # set read orientation for featureCounts (unstranded)
      strand="both"                       # set blast database strandness (unstranded)

    fi

    if [ "$METHOD" != "sensitive" ]; then

      echo "   done... blast unmapped reads on SL sequences..."
      samtools view ${3}_oneEnd_unmapped.bam | awk '{OFS="\t"; print ">"$1"\n"$10}'  > ${3}_oneEnd_unmapped.fasta
      blastn -query ${3}_oneEnd_unmapped.fasta -task blastn -db $SL_db  -outfmt 6 -max_target_seqs 1 -num_threads $threads -word_size 8 -strand $strand > ${3}_blasted.txt 2>${3}_log.txt

      echo "   done... filter significant alignment with qstart==1..."
      awk '$7 == 1 && $10 >= "'"$send"'" && $11 < "0.05" {print $0}' ${3}_blasted.txt | grep "SL1" > ${3}_blasted_SL1.txt
      awk '$7 == 1 && $10 >= "'"$send"'" && $11 < "0.05" {print $0}' ${3}_blasted.txt | grep "SL2" > ${3}_blasted_SL2.txt
      
      echo "   done... trim SL sequences and convert back into fastq..."

      cut -f 1 ${3}_blasted_SL1.txt > ${3}_SL1_IDs.txt
      picard FilterSamReads FILTER=includeReadList COMPRESSION_LEVEL=0 READ_LIST_FILE=${3}_SL1_IDs.txt I=${3}_oneEnd_unmapped.bam O=${3}_SL1.sam 2>> ${3}_log.txt
      samtools sort -n ${3}_SL1.sam > ${3}_SL1_Nsorted.sam
      paste <(samtools view ${3}_SL1_Nsorted.sam | cut -f 1-11) <(cat ${3}_blasted_SL1.txt) | awk '{OFS="\t"; print $1,"4",$3,$4,$5,$6,$7,$8,$9, substr($10, $19+1),  substr($11, $19+1)}' | picard SamToFastq VALIDATION_STRINGENCY=SILENT QUIET=TRUE I=/dev/stdin FASTQ=${3}_SL1_trimmed.fq

      cut -f 1 ${3}_blasted_SL2.txt > ${3}_SL2_IDs.txt
      picard FilterSamReads FILTER=includeReadList COMPRESSION_LEVEL=0 READ_LIST_FILE=${3}_SL2_IDs.txt I=${3}_oneEnd_unmapped.bam O=${3}_SL2.sam 2>> ${3}_log.txt
      samtools sort -n ${3}_SL2.sam > ${3}_SL2_Nsorted.sam
      paste <(samtools view ${3}_SL2_Nsorted.sam | cut -f 1-11) <(cat ${3}_blasted_SL2.txt) | awk '{OFS="\t"; print $1,"4",$3,$4,$5,$6,$7,$8,$9, substr($10, $19+1),  substr($11, $19+1)}' | picard SamToFastq VALIDATION_STRINGENCY=SILENT QUIET=TRUE I=/dev/stdin FASTQ=${3}_SL2_trimmed.fq


    else

      echo "   done... find & cut SL sequences..."

      samtools fixmate ${3}_oneEnd_unmapped.bam - | picard SamToFastq VALIDATION_STRINGENCY=SILENT QUIET=TRUE I=/dev/stdin FASTQ=${3}_oneEnd_unmapped.fq
      cutadapt -g file:$SL_db -O $align_length -m 15 -o ${3}{name}.fq --discard-untrimmed ${3}_oneEnd_unmapped.fq 2>> ${3}_log.txt

      cat ${3}*_SL2_splice_leader*.fq | paste - - - - | sort -n -k2,1 -t. | tr "\t" "\n" > ${3}_SL2_trimmed.fq
      cat ${3}*_SL1_splice_leader*.fq | paste - - - - | sort -n -k2,1 -t. | tr "\t" "\n" > ${3}_SL1_trimmed.fq
      rm ${3}*_SL2_splice_leader*.fq
      rm ${3}*_SL1_splice_leader*.fq

      awk 'NR%4==1 {print substr($1,2); }' ${3}_SL1_trimmed.fq > ${3}_SL1_IDs.txt
      awk 'NR%4==1 {print substr($1,2); }' ${3}_SL2_trimmed.fq > ${3}_SL2_IDs.txt

    fi
    echo "   done... retrieve mapped mates..."

    echo "      [1/2] of SL1-containing reads"
    picard FilterSamReads FILTER=includeReadList READ_LIST_FILE=${3}_SL1_IDs.txt I=${3}_oneEndMapped.bam O=${3}_SL1_mates.bam 2>> ${3}_log.txt
    samtools sort -n ${3}_SL1_mates.bam > ${3}_SL1_mates_Nsorted.bam 
    bedtools bamtofastq -i ${3}_SL1_mates_Nsorted.bam -fq ${3}_SL1_mates.fq

    echo "      [2/2] of SL2-containing reads"
    picard FilterSamReads FILTER=includeReadList READ_LIST_FILE=${3}_SL2_IDs.txt I=${3}_oneEndMapped.bam O=${3}_SL2_mates.bam 2>> ${3}_log.txt
    samtools sort -n ${3}_SL2_mates.bam > ${3}_SL2_mates_Nsorted.bam 
    bedtools bamtofastq -i ${3}_SL2_mates_Nsorted.bam -fq ${3}_SL2_mates.fq

    echo "   done... remap with hisat2"
    hisat2 -p $threads --no-discordant --no-softclip --min-intronlen 20 --max-intronlen 5000 --rna-strandness $paired_orientation -x $index -1 ${3}_SL1_mates.fq -2 ${3}_SL1_trimmed.fq | samtools view -b -F 260 -f 2 > ${3}_SL1_remapped.bam
    hisat2 -p $threads --no-discordant --no-softclip --min-intronlen 20 --max-intronlen 5000 --rna-strandness $paired_orientation -x $index -1 ${3}_SL2_mates.fq -2 ${3}_SL2_trimmed.fq | samtools view -b -F 260 -f 2 > ${3}_SL2_remapped.bam

    echo "   done... summarize..."
    featureCounts -s -p $featureCounts_S -g gene_id -T 4 -a $gene_annotation -o ${3}_counts.txt ${3}_SL1_remapped.bam ${3}_SL2_remapped.bam 2>> ${3}_log.txt

  fi


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# SINGLE-END MODE
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

else
  if [ $# != 2 ]; then
     echo "${usage}"
     exit 1

  else
    echo ""
    echo "SL-quant.sh: SINGLE-END MODE [$METHOD]"
    echo "   [unmapped.bam] = $1"
    echo "   [ouput_dir/base] = $2"
    echo "   [SL blast database] = $SL_db"
    echo "   [gene annotation] = $gene_annotation"
    echo "   [read orientation] = $single_orientation";echo ""

 
    if [ "$single_orientation" == "R" ]; then
      strand="plus"                                 # set blast database strandness (normal)
      featureCounts_S=1                             # set featureCounts strandness (stranded)    

    elif [ "$single_orientation" == "F" ]; then
      strand="minus"                                # set blast database strandness (reverse)
      featureCounts_S=1                             # set featureCounts strandness (stranded) 

    else
      strand="both"                                 # set blast database strandness (unstranded)
      featureCounts_S=0                             # set featureCounts strandness (unstranded)
    fi
    
    paired_entries=$(samtools view -c -f 1 $1)
    
    if [ $paired_entries -gt 0 ]; then
      echo ""; echo "WARNING : there are ${paired_entries} 'paired in sequencing' reads in ${1}. Consider running the script in paired-end mode. Only R2 reads will be considered in single-end mode."; echo ""
      samtools view -b -f 128 ${1} > ${1}_R2.bam
      input="${1}_R2.bam"
    
    else
      input="${1}"
    fi


    if [ "$METHOD" != "sensitive" ]; then

      echo "   blast unmapped reads on SL sequences..."
      samtools view -f 4 $input | awk '{OFS="\t"; print ">"$1"\n"$10}' | blastn -db $SL_db -outfmt 6 -max_target_seqs 1 -num_threads 4 -word_size 8 -strand $strand > ${2}_blasted.txt 2>${2}_log.txt
    

      echo "   done... filter..."
      awk '$7 == 1 && $10 >= "'"$send"'" && $11 < "0.05" {print $0}' ${2}_blasted.txt | grep "SL1" > ${2}_blasted_SL1.txt
      awk '$7 == 1 && $10 >= "'"$send"'" && $11 < "0.05" {print $0}' ${2}_blasted.txt | grep "SL2" > ${2}_blasted_SL2.txt

      echo "   done... retrieve SL-containing reads"

      echo "      [1/2] of SL1-containing reads"
      cut -f 1 ${2}_blasted_SL1.txt > ${2}_SL1_IDs.txt
      picard FilterSamReads FILTER=includeReadList VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 QUIET=TRUE READ_LIST_FILE=${2}_SL1_IDs.txt I=$input O=${2}_SL1.sam


      echo "      [2/2] of SL2-containing reads"
      cut -f 1 ${2}_blasted_SL2.txt > ${2}_SL2_IDs.txt
      picard FilterSamReads FILTER=includeReadList VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 QUIET=TRUE READ_LIST_FILE=${2}_SL2_IDs.txt I=$input O=${2}_SL2.sam

      echo "   done... trim SL-containing reads and convert to fastq"
      paste <(samtools view ${2}_SL1.sam | cut -f 1-11) <(cat ${2}_blasted_SL1.txt) | awk '{OFS="\t"; print $1,"4",$3,$4,$5,$6,$7,$8,$9, substr($10, $19+1),  substr($11, $19+1)}' | picard SamToFastq VALIDATION_STRINGENCY=SILENT QUIET=TRUE I=/dev/stdin FASTQ=${2}_SL1_trimmed.fq
      paste <(samtools view ${2}_SL2.sam | cut -f 1-11) <(cat ${2}_blasted_SL2.txt) | awk '{OFS="\t"; print $1,"4",$3,$4,$5,$6,$7,$8,$9, substr($10, $19+1),  substr($11, $19+1)}' | picard SamToFastq VALIDATION_STRINGENCY=SILENT QUIET=TRUE I=/dev/stdin FASTQ=${2}_SL2_trimmed.fq


    else

      SL_db="data/blast_db/SL_RV.fa"         # path to SL sequence database (for blast).

      echo "   find & cut SL sequences..."

      samtools fixmate ${input} - | picard SamToFastq VALIDATION_STRINGENCY=SILENT QUIET=TRUE I=/dev/stdin FASTQ=${2}.fq
      cutadapt -g file:$SL_db -O $align_length -m 15 -o ${2}{name}.fq --discard-untrimmed ${2}.fq 2>> ${2}_log.txt

      cat ${2}*_SL2_splice_leader*.fq | paste - - - - | sort -n -k2,1 -t. | tr "\t" "\n" > ${2}_SL2_trimmed.fq
      cat ${2}*_SL1_splice_leader*.fq | paste - - - - | sort -n -k2,1 -t. | tr "\t" "\n" > ${2}_SL1_trimmed.fq
      rm ${2}*_SL2_splice_leader*.fq
      rm ${2}*_SL1_splice_leader*.fq

    fi
 
    echo "   done... map with hisat2"
    hisat2 -p $threads --no-discordant --no-softclip --min-intronlen 20 --max-intronlen 5000 --rna-strandness $single_orientation -x $index -U ${2}_SL1_trimmed.fq | samtools view -b -F 260 > ${2}_SL1_remapped.bam
    hisat2 -p $threads --no-discordant --no-softclip --min-intronlen 20 --max-intronlen 5000 --rna-strandness $single_orientation -x $index -U ${2}_SL2_trimmed.fq | samtools view -b -F 260 > ${2}_SL2_remapped.bam

    echo "   done... summarize..."
    featureCounts -s 1 -g gene_id -T 4 -a $gene_annotation -o ${2}_counts.txt ${2}_SL1_remapped.bam ${2}_SL2_remapped.bam 2>> ${2}_log.txt
  fi

fi
echo "   All done ! Have a nice day :)"
echo ""
