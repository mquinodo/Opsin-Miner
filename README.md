# Opsin-Miner
Tool to find coding variants, deletions and pathogenic haplotypes in OPN1LW, OPN1MW and OPN1MW2 using exome sequencing data

This software was written by Mathieu Quinodoz in the group of Prof. Rivolta from the IOB in Basel, Switzerland. It is presented at XYZ. It was developped on Ubuntu 20.04.5 LTS (GNU/Linux 5.4.0-131-generic x86_64).

## Prerequisites (version recommended)
+ BCFTools [[Link](https://samtools.github.io/bcftools/howtos/install.html)] (>= v1.9)
+ BEDTools [[Link](https://bedtools.readthedocs.io/en/latest/content/installation.html)] (>= v2.25.0)
+ R [[Link](https://cran.r-project.org/mirrors.html)] (>= v3.2.0)
+ BWA [[Link](https://sourceforge.net/projects/bio-bwa/files)] (>=0.7.17)
+ Java (in order to run Picard and GATK)
+ Picard jar file [[Link](https://github.com/broadinstitute/picard)] (>= v2.23)
+ GATK jar file [[Link](https://github.com/broadinstitute/gatk/releases)] (>= v4.0.0.0)

## Installation
The tool does not require compilation. It will download reference genome files from UCSC during the first run.

## Usage
This tool is made to be run by batch of samples ideally sequenced in the same batch or at least with same capture kit and sequencer. This is required for proper comparison between samples to detect deletions.

The main script 02_opsins-main.sh takes as input a tsv file with 2 or 3 columns: ID, forward FASTQ file and reverse FASTQ file (for paired-end). From these FASTQ files, Opsin-Miner will determine possible deletions, detect pathogenic haplotypes and coding variants.
It is called with bash:
```
bash 02_opsins-main.sh --script script_folder --out output_folder --input input_file.tsv --picard picard_jar_location --gatk gatk_jar_location
```
The approximate computation time per sample is few minutes per sample.

#### Required arguments
Option  | Value  | Description 
--- | --- | ---
--script | STRING | Folder where the scripts are (e.g. XYZ/Opsin-Miner-main)
--out | STRING | Folder for outputs
--input | STRING | Tab-delimited text file with 2 or 3 columns: ID, forward FASTQ file and eventually reverse FASTQ file
--picard | STRING | JAR file to run Picard
--gatk | STRING | JAR file to run GATK

