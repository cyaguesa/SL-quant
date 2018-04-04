#!/bin/bash

chmod +x SL-quant.sh map_reads_modENCODE.sh map_reads.sh SL_sites.sh SL_cutadapt.sh
mkdir data/ce10_hisat2_index
mkdir SL-quant_results
echo ""

echo 'download C. elegans annotation...'
curl ftp://ftp.wormbase.org/pub/wormbase/releases/WS262/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS262.canonical_geneset.gtf.gz -o data/genes.gtf.gz
gunzip data/genes.gtf.gz
echo ""

echo 'download C. elegans genome...'
curl ftp://ftp.wormbase.org/pub/wormbase/releases/WS262/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS262.genomic.fa.gz -o data/ce10_hisat2_index/genome.fa.gz
gunzip data/ce10_hisat2_index/genome.fa.gz
echo ""

echo 'build hisat2 index...'
hisat2-build -p 4 data/ce10_hisat2_index/genome.fa data/ce10_hisat2_index/genome
echo ""

echo 'build samtools index...'
samtools faidx data/ce10_hisat2_index/genome.fa
echo ""

echo 'build blast database...'
makeblastdb -dbtype nucl -in data/blast_db/SL.fasta
echo ""

echo 'done ! You can now test SL-quant !'
