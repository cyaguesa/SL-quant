#!/bin/bash

# Implementation of the SL-containing reads identification strategy used in (Tourasse et al., 2017) [doi:10.1101/gr.224626.117]

# REQUIREMENTS (see installation instruction if unmet):      VERSION USED       NOTE

# - picard installed and in your path                        2.9.0
# - hisat2 installed and in your path                        2.0.5             
# - cutadapt installed and in your pat                       1.14              
# - samtools installed and in your path                      1.5

# by Carlo Yague-Sanz, 2017


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# PARAMETERS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

set -e                                 # stops script at the first error.
SL_db="data/blast_db/SL_RV.fa"         # path to SL sequence (for cutadapt).
#SL_db="data/blast_db/SL.fasta"         # path to SL sequence database (for blast).
index="data/ce10_hisat2_index/genome"  # genome index file.
align_length=5                         # length of the cutadapt alignment (default = 5).
threads=4                              # number of threads to use.
read_orientation="R"                   # read orientation (default="R", reversely stranded).                     
                     

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# "CUTADAPT" MODE
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  if [ $# != 2 ]; then
    echo ""; echo "   Usage : ./SL_cutadapt.sh [unmapped.bam] [ouput_dir/base]"; echo ""
    echo "   Example : ./SL_cutadapt.sh data/reads/modENCODE/modENCODE_4594/unmapped.bam SL-quant_results/modENCODE_4594_cutadapt/"
    echo "   Example : ./SL_cutadapt.sh data/reads/SRR1585277/SRR1585277/unmapped.bam_R2.bam SL-quant_results/SRR1585277_single_cutadapt/"

    exit 1

  else
    echo ""
    echo "SL-cutadapt.sh: 'CUTADAPT' MODE"
    echo "   [unmapped.bam] = $1"
    echo "   [ouput_dir/base] = $2"
    echo "   [SL fasta] = $SL_db"
    echo "   [read orientation] = $read_orientation"
    echo "   [alignment min length] = $align_length";echo ""

    echo "   convert unmapped reads to fastq ..."
    mkdir -p "$(dirname ${2}_)"
    bedtools bamtofastq -i ${1} -fq ${2}_unmapped.fq

    echo "   done... trim SL-sequences..."

    cutadapt -g file:$SL_db -O $align_length -m 15 -o ${2}{name}.fq --discard-untrimmed ${2}_unmapped.fq

    cat ${2}*_SL2_splice_leader*.fq > ${2}_SL2_merged.fq
    cat ${2}*_SL1_splice_leader*.fq > ${2}_SL1_merged.fq
    rm ${2}*_SL2_splice_leader*.fq
    rm ${2}*_SL1_splice_leader*.fq

    echo "   done... remap with hisat2"
    hisat2 -p $threads --no-discordant --no-softclip --min-intronlen 20 --max-intronlen 5000 --rna-strandness $read_orientation -x $index -U ${2}_SL2_merged.fq | samtools view -b -F 260 > ${2}_SL2_remapped.bam
    hisat2 -p $threads --no-discordant --no-softclip --min-intronlen 20 --max-intronlen 5000 --rna-strandness $read_orientation -x $index -U ${2}_SL1_merged.fq | samtools view -b -F 260 > ${2}_SL1_remapped.bam

   fi

echo "   All done ! Have a nice day :)"
echo ""


