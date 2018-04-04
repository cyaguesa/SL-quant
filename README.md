# SL-quant

## A pipeline to quantify SL trans-splicing events from RNA-seq data

_SL-quant_ is a bash pipeline that adapts to paired-end and single-end RNA-seq data to accurately quantify splice-leader (SL) trans-splicing events by genes in the nematode _C. elegans_. It is designed to work downstream of read mapping and takes the reads left unmapped as primary input. _SL-quant_ completes under 15 minutes on a basic desktop computer for typical RNA-seq libraries.

Detailed description and validation of the pipeline are reported in the manuscript [still under consideration for publication; reference will be added once published].

For support, questions or requests, please contact: carlo.yaguesanz@unamur.be 

- [Installation & quick start](https://github.com/cyaguesa/SL-quant/blob/master/README.md#installation--quick-start)
- [Detailed usage](https://github.com/cyaguesa/SL-quant/blob/master/README.md#detailed-usage)
- [Adaptation to other species](https://github.com/cyaguesa/SL-quant/blob/master/README.md#adaptation-to-other-species)
- [Identification of trans-splice sites](https://github.com/cyaguesa/SL-quant/blob/master/README.md#identification-of-trans-splice-sites)
- [Reproduce the analysis from the manuscript](https://github.com/cyaguesa/SL-quant/blob/master/README.md#reproduce-the-analysis-from-the-manuscript)

## Installation & quick start

_SL-quant_ comes as a simple bash script that works on macOS and Linux systems. However, the following dependencies need to be installed and set in your PATH:

- [blastn](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) from the blast+ suite (2.6.0 or newer)
- [samtools](http://samtools.sourceforge.net/) (1.5 or newer)
- [picard-tools](http://broadinstitute.github.io/picard/) (2.9.0 or newer)
- [featureCounts](http://subread.sourceforge.net/) from the subread package. (1.5.0 or newer)
- [bedtools](http://bedtools.readthedocs.io/en/latest/) (2.26.0 or newer)
- [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html) (1.14 or newer)
- [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) (2.0.5 or newer)


There are many ways to install these but the easiest might be through a package manager such as brew. 

#### Easy dependencies installation with homebrew

###### Install homebrew on macOS
Follow instructions on the [brew homepage](https://brew.sh/) or simply paste this command in the Terminal:

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

###### Install linuxbrew on linux systems
Follow instructions on the [linuxbrew homepage](https://linuxbrew.sh/).

###### Install dependencies with brew (macOS and linux)
Paste the following commands in the Terminal:

    brew update
    brew install blast
    brew install samtools
    brew install picard-tools
    brew install bedtools
    brew install hisat2
    brew install cutadapt
    brew tap jrderuiter/tap
    brew install subread

#### Setting-up and testing SL-quant

###### Clone SL-quant
To install _SL-quant_, clone this repository using git:

    git clone https://github.com/cyaguesa/SL-quant

or alternatively, use the green downoad button at the top of the page.

###### Set-up SL-quant
To complete SL-quant set-up, use this script to create directories, make the other scripts executable, download the _C. elegans_ genome sequence and build the indexes. This can take a few minutes.

    cd SL-quant  # or 'cd ~/Desktop/SL-quant' if it is cloned on your desktop
    chmod +x set-up.sh
    ./set-up.sh

###### Test SL-quant
We provide a small test dataset (1000 reads) for testing SL-quant in a few seconds.

    time ./SL-quant.sh data/reads/test_unmapped.bam SL-quant_results/test
    real	0m2.718s
    user	0m3.235s
    sys	0m0.589s
    
###### Basic usage

    ./SL-quant.sh -h

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


## Detailed usage

#### A word about the SL-quant input
To minimize _SL-quant_ running time, its primary input is limited to reads susceptible to originate from trans-spliced RNA fragments, that is, unmapped reads in bam format. This implies that a first round of mapping to the _C. elegans_ genome or transcriptome must precede the use of _SL-quant_. It must be performed end-to-end (without soft-clipping) in order to make sure that reads originating from trans-spliced RNA fragments do not map. This is the default behaviour of many mappers (bowtie2, tophat2, BBMap, …) but for others, such as STAR or HiSAT2, soft-clipping should be disabled. Beside this specification, any coordinate-sorted bam file containing unmapped reads can be fed into _SL-quant_, making it particularly well suited for posterior analysis of old data.

#### Basic paramaters

###### paired-end mode (-p --paired)
In the case paired-end data is available, we provide an optimized paired-end mode that requires an additional parameter -m --mapped, referring to the mapped reads of the library. These mapped reads will be used to prefilter the unmapped reads and to refine the mapping after the SL-sequence trimming.

Example: `./SL-quant.sh -p -m data/reads/test_mapped.bam data/reads/test_unmapped.bam SL-quant_results/test`

###### sensitive mode (-s --sensitive)
In sensitive mode, the criteria to identify SL sequences are less stringent, increasing sensitivity (at the cost of some specificity). It can be used in combinaison with the paired-end mode or in single-end mode.

Example 1 : `./SL-quant.sh -s data/reads/test_unmapped.bam SL-quant_results/test`

Example 2 : `./SL-quant.sh -p -s -m data/reads/test_mapped.bam data/reads/test_unmapped.bam SL-quant_results/test`


#### Advanced paramaters

At the beginning of the SL-quant.sh file, there are a few additional parameters that can be modified for advanced users.

- `set -e` This means that the script will stop at the first error (which is usually for the best).
- `SL_db="data/blast_db/SL.fasta"` The path to SL sequence database for blast.
- `gene_annotation="data/genes.gtf"`The annotation file for the summarization step.
- `index="data/ce10_hisat2_index/genome"`The genome index for hisat2.
- `paired_orientation="FR"`The read orientation for paired-end mode (ignored in single-end mode). Value={"FR" (default), "RF", "unstranded"}
- `single_orientation="R"`The read orientation for single-end mode (ignored in paired-end mode). Value={"F" (stranded), "R" (reversely stranded), "unstranded"}
- `threads=4` The number of threads to use.
- `send=22` The minimal end of alignment in 'subject' threshold for blast (this should be equal to the length of the SL sequence).
- `align_length=5` The minimal length of the cutadapt alignment (default = 5).

#### Output files

###### Counts results
The final output file named `[output_dir/base]_counts.txt` is a tab delimited file generated by featureCounts containing the number of SL1 and SL2 trans-splicing events by genes. A summary of the quantification is available in the `[output_dir/base]_counts.txt.summary` file.

###### Intermediate files
Those files are generated during the _SL-quant_ process. When the file name ends by SL1, there is a similar file ending with SL2 for the SL2-containing reads.
- `[output_dir/base]_blasted.txt` countains the raw results of the blast of the unmapped reads to the SL sequences.
- `[output_dir/base]_blasted_SL1.txt` countains blast results for the SL1-containing reads.
- `[output_dir/base]_SL1_IDs.txt` countains the read IDs of the SL1-containing reads.
- `[output_dir/base]_SL1.sam` countains the SL1-containing reads in SAM format (single-end mode only).
- `[output_dir/base]_SL1_trimmed.fq` countains the trimmed SL1-containing reads in fastq format (single-end mode only).
- `[output_dir/base]_SL1_remapped.bam` countains the trimmed SL1-containing reads mapped on the genome in bam format (single-end mode only).
- `[output_dir/base]_oneEnd_unmapped.fasta` countains unmapped reads (after prefiltering) in fasta format (paired-end mode only).
- `[output_dir/base]_oneEndMapped.bam` countains the mates of the unmapped reads (after prefiltering) in bam format (paired-end mode only).

###### log file
`[output_dir/base]_log.txt` countains various warning/outputs generated during the _SL-quant_ process.

## Adaptation to other species
While SL-quant was developed for and tested on _C.elegans_ data, many other species do SL-trans-splicing. The analysis of trans-splicing events in such species is possible with SL-quant with the following adaptations. Don't hesitate to contact us if you need any help implementing those changes.

#### Change the SL sequences in the blast database
1- Find the SL sequences (not the full SL RNA sequences, only the part that will be trans-spliced to the mRNA) and save it in a fasta file `SL_my_species.fasta` in the `data/blast_db` directory. You should include the characters "SL1" and/or "SL2" in the header of the fasta sequence.

2- Build the new blast database:

    makeblastdb -dbtype nucl -in data/blast_db/SL_my_species

3- Replace the value of the `SL_db` parameter in the `SL-quant.sh` script by `"data/blast_db/SL_my_species.fasta"`.

4- Replace the value of the `send` parameter in the `SL-quant.sh` script by the length of the new SL sequence(s).

#### Change the reference genome index.
1- Download the reference genome for your species of interest and save it as a fasta file `genome_my_species.fa` in a new `data/index_my_species` directory.

2- Build the hisat2 index:

    hisat2-build data/index_my_species/genome_my_species.fa data/index_my_species/genome_my_species

3- Replace the value of the `index` parameter in the `SL-quant.sh` script by `"data/index_my_species/genome_my_species.fa"`.

#### Change the gene annotation file.
1 - For now, SL-quant only supports [SAF annotation files](http://bioinf.wehi.edu.au/featureCounts/). Support for .gtf files is planned in the near future. Download or create one of those file `genes_my_species.SAF` for your species and save it into the `data` directory.

2- Replace the value of the `gene_annotation` parameter in the `SL-quant.sh` script by `data/genes_my_species.SAF`.

## Identification of trans-splice sites

We designed SL-quant with the idea of quantifying SL trans-splicing events by genes but it is also possible to identify trans-splice sites at single nucleotide resolution from the output. Indeed, in single-end mode, the 5' end of the reads mapped after SL sequence trimming correspond to the position of the trans-splice sites. The folowing lines describe such analysis applied to the SL1 trans-splicing only.

#### sort remapped bam files

    samtools sort SL-quant_results/test_SL1_remapped.bam -o SL-quant_results/test_SL1_remapped_sorted.bam

#### get 5' end positions of reads (strand-specific)

    bedtools genomecov -ibam SL-quant_results/test_SL1_remapped_sorted.bam -dz -5 -strand + > SL-quant_results/test_TS_SL1_sites_fwd.tab
    bedtools genomecov -ibam SL-quant_results/test_SL1_remapped_sorted.bam -dz -5 -strand - > SL-quant_results/test_TS_SL1_sites_fwd.tab

#### see transpliced site

    head -n 1 SL-quant_results/test_TS_SL1_sites_fwd.tab
    chrI	6789739     1   # 1 trans-splicing event at position 6789739 of strand '+' of chrI

## Reproduce the analysis from the manuscript.
To reproduce the full analysis presented in our manuscript from the raw data, [R](https://www.r-project.org/) and [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) should be installed. Trimmomatic can be installed with brew:

    brew tap brewsci/bio
    brew install trimmomatic

### Get and map data
We provide two bash scripts to map the reads for the paired-end (`map_reads.sh`) and the single-end dataset (`./map_reads_modENCODE.sh`).

###### paired-end dataset (SRR1585277)
    mkdir data/reads/SRR1585277
    curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR158/007/SRR1585277/SRR1585277_1.fastq.gz -o data/reads/SRR1585277/SRR1585277_1.fastq.gz
    curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR158/007/SRR1585277/SRR1585277_2.fastq.gz -o data/reads/SRR1585277/SRR1585277_2.fastq.gz
    ./map_reads.sh

###### single-end dataset (modENCODE_4594)
    mkdir data/reads/modENCODE
    curl ftp://data.modencode.org/all_files/cele-raw-1/4594_SRR125481.fastq.gz -o data/reads/modENCODE/modENCODE_4594.fastq.gz
    ./map_reads_modENCODE.sh data/reads/modENCODE/modENCODE_4594.fastq.gz

###### generate random reads
    bedtools random -l 50 -seed 0 -n 1000000 -g data/chrom_summary.txt > data/reads/random.bed
    bedtools getfasta -fi data/ce10_hisat2_index/genome.fa -bed data/reads/random.bed > data/reads/random.fa

### Run SL-quant

###### paired-end dataset (SRR1585277)
In paired-end mode (set SINGLE parameter to "paired")

    ./SL-quant.sh data/reads/SRR1585277/SRR1585277/accepted_hits_sorted.bam data/reads/SRR1585277/SRR1585277/unmapped.bam SL-quant_results/SRR1585277_paired"

In single-end mode (set SINGLE parameter to "single")

    ./SL-quant.sh data/reads/SRR1585277/SRR1585277/unmapped.bam SL-quant_results/SRR1585277_single"

###### single-end dataset (modENCODE_4594)
    ./SL-quant.sh data/reads/modENCODE/modENCODE_4594/unmapped.bam SL-quant_results/modENCODE_4594

###### blast random reads
    blastn -query data/reads/random.fa -db data/blast_db/SL.fasta -outfmt 6 -max_target_seqs 1 -num_threads 4 -word_size 8 > SL-quant_results/random_blasted.txt

## Analyse the data and make the plots
An R script (`analyse_SL.R`) is provided to reproduce the analysis presented in the manuscript.
