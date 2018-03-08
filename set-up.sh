#!/bin/bash

chmod +x SL-quant.sh map_reads_modENCODE.sh map_reads.sh
mkdir data/ce10_bowtie2_index
mkdir SL-quant_results
echo 'download C. elegans genome...'
curl ftp://ftp.wormbase.org/pub/wormbase/releases/WS220/species/c_elegans/c_elegans.WS220.genomic.fa.gz -o data/ce10_bowtie2_index/genome.fa.gz
echo 'rename chromosomes...'
gunzip data/ce10_bowtie2_index/genome.fa.gz
sed 's/CHROMOSOME_/chr/g' data/ce10_bowtie2_index/genome.fa | sed 's/chrMtDNA/chrM/g' > data/ce10_bowtie2_index/temp.fa
mv data/ce10_bowtie2_index/temp.fa data/ce10_bowtie2_index/genome.fa
rm data/ce10_bowtie2_index/temp.fa
echo 'build bowtie2 index...'
bowtie2-build --threads 4 data/ce10_bowtie2_index/genome.fa data/ce10_bowtie2_index/genome
echo 'build hisat2 index...'
hisat2-build
echo 'build samtools index...'
samtools faidx data/ce10_bowtie2_index/genome.fa
echo 'done ! You can now test SL-quant !'
