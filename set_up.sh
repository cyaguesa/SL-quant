#!/bin/bash

echo 'download C. elegans genome...'
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS220/species/c_elegans/c_elegans.WS220.genomic.fa.gz -P data/ce10_bowtie2_index
echo 'gunzip it...'
gunzip data/ce10_bowtie2_index/c_elegans.WS220.genomic.fa.gz
mv data/ce10_bowtie2_index/c_elegans.WS220.genomic.fa data/ce10_bowtie2_index/genome.fa
echo 'build bowtie2 index...'
bowtie2-build data/ce10_bowtie2_index/genome.fa data/ce10_bowtie2_index/genome
echo 'build samtools index...'
samtools faidx data/ce10_bowtie2_index/genome.fa
echo 'done ! You can now test SL-quant !'