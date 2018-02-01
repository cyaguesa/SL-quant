#!/bin/bash

mkdir data/ce10_bowtie2_index
mkdir SL-quant_results
echo 'download C. elegans genome...'
curl ftp://ftp.wormbase.org/pub/wormbase/releases/WS220/species/c_elegans/c_elegans.WS220.genomic.fa.gz -o data/ce10_bowtie2_index/genome.fa.gz
echo 'gunzip it...'
gunzip data/ce10_bowtie2_index/genome.fa.gz
echo 'build bowtie2 index...'
bowtie2-build data/ce10_bowtie2_index/genome.fa data/ce10_bowtie2_index/genome
echo 'build samtools index...'
samtools faidx data/ce10_bowtie2_index/genome.fa
echo 'done ! You can now test SL-quant !'
