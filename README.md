# SL-quant

## A pipeline to quantify SL trans-splicing events from RNA-seq data

_SL-quant_ is a bash pipeline that adapts to paired-end and single-end RNA-seq data to accurately quantify splice-leader (SL) trans-splicing events by genes in the nematode _C. elegans_. It is designed to work downstream of read mapping and takes the reads left unmapped as primary input. _SL-quant_ completes under 15 minutes on a basic desktop computer for typical RNA-seq libraries.

Detailed description and validation of the pipeline are reported in the manuscript [still under consideration for publication; reference will be added once published].

For support, questions or requests, please contact: carlo.yague-sanz@unamur.be 

## Installation & quick start

_SL-quant_ comes as a simple bash script that works on macOS and Linux systems. However, the following dependencies need to be installed and set in your PATH:

- [blastn](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) from the blast+ suite
- [samtools](http://samtools.sourceforge.net/)
- [picard-tools](http://broadinstitute.github.io/picard/)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [featureCounts](http://subread.sourceforge.net/) from the subread package.

There are many ways to install these but the easiest might be through a package manager such as brew. 

#### Easy dependencies installation on macOS with homebrew

###### Install homebrew
Follow instructions on the [brew homepage](https://brew.sh/) or simply paste this command in the Terminal:

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

###### Install dependencies with brew
Paste the following commands in the Terminal:

    brew update
    brew install blast
    brew install samtools
    brew install picard-tools
    brew install bowtie2
    brew tap jrderuiter/tap
    brew install subread

#### Easy dependencies installation on Linux with linuxbrew

###### Install linuxbrew
Follow instructions on the [linuxbrew homepage](https://linuxbrew.sh/).

###### Install dependencies with linuxbrew
Paste the following commands in the Terminal:

    brew update
    brew install blast
    brew install samtools
    brew install picard-tools
    brew install bowtie2
    brew tap jrderuiter/tap
    brew install subread

#### Setting-up and testing SL-quant

###### Clone SL-quant
To install _SL-quant_, clone this repository using git:

    git clone https://github.com/cyaguesa/SL-quant

or alternatively, use the green downoad button at the top of the page.

###### Set-up SL-quant
To complete SL-quant set-up, use this script to download the _C. elegans_ genome sequence and build the indexes:

    cd SL-quant
    chmod +x SL-quant.sh set-up.sh
    ./set-up.sh

###### Test SL-quant
We provide a small test dataset (1000 reads) for testing SL-quant with its default parameters.

    ./SL-quant.sh data/reads/test_unmapped.bam SL-quant_results/test

## Detailed usage

#### A word about the SL-quant input
To minimize _SL-quant_ running time, its primary input is limited to reads susceptible to originate from trans-spliced RNA fragments, that is, unmapped reads in bam format. This implies that a first round of mapping to the _C. elegans_ genome or transcriptome must precede the use of _SL-quant_. It must be performed end-to-end (without soft-clipping) in order to make sure that reads originating from trans-spliced RNA fragments do not map. This is the default behaviour of many mappers (bowtie2, tophat2, BBMap, …) but for others, such as STAR or HiSAT2, soft-clipping should be disabled. Beside this specification, any coordinate-sorted bam file containing unmapped reads can be fed into _SL-quant_, making it particularly well suited for posterior analysis of old data.

#### Basic paramaters

###### single-end mode (default)
In its default mode _SL-quant_ takes two parameters, the unmapped reads in bam format and the base name (path/name) for the outputs files:

Usage: `./SL-quant.sh [unmapped.bam] [ouput_dir/base]`

Example: `./SL-quant.sh data/reads/test_unmapped.bam SL-quant_results/test`

###### paired-end mode
In the case paired-end data is available, we provide an optimized paired-end mode (see advanced parameters) that require three parameters:

Usage: `./SL-quant.sh [mapped.bam] [unmapped.bam] [ouput_dir/base]`

Example: `./SL-quant.sh data/reads/test_mapped.bam data/reads/test_unmapped.bam SL-quant_results/test_paired`

#### Advanced paramaters

At the beginning of the SL-quant.sh file, there are a few additional parameters that can be modified for advanced users.

- `set -e` This means that the script will stop at the first error (which is usually for the best).
- `SINGLE="single"` Set to "single" for single-end mode (default). Any other value triggers paired-end mode.
- `SL_db="data/blast_db/SL.fasta"` The path to SL sequence database for blast.
- `gene_annotation="data/genes.gtf"`The annotation file for the summarization step.
- `index="data/ce10_bowtie2_index/genome"`The genome index for bowtie2 (only required on single-end mode).
- `paired_orientation="fr-firststrand"`The read orientation for paired-end mode (ignored in single-end mode). Should be set to either "fr-firststrand" (default), "fr-secondstrand" or "fr-unstrand".
- `single_orientation="stranded"`The read orientation for single-end mode (ignored in paired-end mode). Should be set to either "stranded" (default), "reversely_stranded" or "unstranded".

#### Output files

###### Counts results
The final ouput file named `[ouput_dir/base]_counts.txt` is a tab delimited file generated by featureCounts containing the number of SL1 and SL2 trans-splicing events by genes. A summary of the quantification is available in the `[ouput_dir/base]_counts.txt.summary`  file.

###### Intermediate files
Those files are generated during the _SL-quant_ process. When the file name ends by SL1, there is a similar file ending with SL2 for the SL2-containing reads.
- `[ouput_dir/base]_blasted.txt` countains the raw results of the blast of the unmapped reads to the SL sequences.
- `[ouput_dir/base]_blasted_SL1.txt` countains blast results for the SL1-containing reads.
- `[ouput_dir/base]_SL1_IDs.txt` countains the read IDs of the SL1-containing reads.
- `[ouput_dir/base]_SL1.sam` countains the SL1-containing reads in SAM format (single-end mode only).
- `[ouput_dir/base]_SL1_trimmed.fq` countains the trimmed SL1-containing reads in fastq format (single-end mode only).
- `[ouput_dir/base]_SL1_remapped.bam` countains the trimmed SL1-containing reads mapped on the genome in bam format (single-end mode only).
- `[ouput_dir/base]_oneEnd_unmapped.fasta` countains unmapped reads (after prefiltering) in fasta format (paired-end mode only).
- `[ouput_dir/base]_oneEndMapped.bam` countains the mates of the unmapped reads (after prefiltering) in bam format (paired-end mode only).

###### log file
`[ouput_dir/base]_log.txt` countains various warning/outputs generated during the _SL-quant_ process.

## Reproduce the analysis of the manuscript.
To reproduce the full analysis from the raw data, [R](https://www.r-project.org/), [bedtools](http://bedtools.readthedocs.io/en/latest/), [tophat2](https://ccb.jhu.edu/software/tophat/index.shtml) and [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) should be installed.

### Get and map data
We provide two bash scripts to map the reads for the paired-end (`map_reads.sh`) and the single-end dataset (`./map_reads_modENCODE.sh`).

###### paired-end dataset (SRR1585277)
    cd ~/Desktop/SL-quant
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR158/007/SRR1585277/SRR1585277_1.fastq.gz -P data/reads/SRR1585277
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR158/007/SRR1585277/SRR1585277_2.fastq.gz -P data/reads/SRR1585277
    ./map_reads.sh

###### single-end dataset (modENCODE_4594)
    wget ftp://data.modencode.org/all_files/cele-raw-1/4594_SRR125481.fastq.gz -P data/reads/modENCODE
    ./map_reads_modENCODE.sh

###### generate random reads
    bowtie2-inspect  data/ce10_bowtie2_index/genome > data/ce10_bowtie2_index/genome.fa
    bedtools random -l 50 -seed 0 -n 1000003 -g data/chrom_summary.txt > data/reads/random.bed
    bedtools getfasta -fi data/ce10_bowtie2_index/genome.fa -bed data/reads/random.bed > data/reads/random.fa

### Run SL-quant

###### paired-end dataset (SRR1585277)
In paired-end mode (set SINGLE parameter to "paired")

    ./SL-quant.sh data/reads/SRR1585277/accepted_hits_sorted.bam data/reads/SRR1585277/unmapped.bam SL-quant_results/SRR1585277_paired"

In single-end mode (set SINGLE parameter to "single")

    ./SL-quant.sh data/reads/SRR1585277/unmapped.bam SL-quant_results/SRR1585277_single"

###### single-end dataset (modENCODE_4594)
    ./SL-quant.sh data/reads/modENCODE/modENCODE_4594/unmapped.bam SL-quant_results/modENCODE_4594

###### blast random reads
    blastn -query data/reads/random.fa -db data/blast_db/SL.fasta -outfmt 6 -max_target_seqs 1 -num_threads 4 -word_size 8 > SL-quant_results/random_blasted.txt

## Analyse the data and make the plots
An R script (`analyse_SL.R`) is provided to reproduce the analysis presented in the analysis.
