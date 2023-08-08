# Opsin-Miner
Tool to find coding variants, deletions and pathogenic haplotypes in OPN1LW, OPN1MW and OPN1MW2 using exome sequencing data

This software was written by Mathieu Quinodoz in the group of Prof. Rivolta from the IOB in Basel, Switzerland in collaboration with Dr. Susanne Kohl and Prof. Bernd Wissinger from University of Tübingen, Germany.

It was developped on Ubuntu 20.04.5 LTS (GNU/Linux 5.4.0-131-generic x86_64).

## Prerequisites (version recommended)
+ SAMtools [[Link](http://www.htslib.org/download/)] (>= v1.10)
+ BCFTools [[Link](https://samtools.github.io/bcftools/howtos/install.html)] (>= v1.9)
+ BEDTools [[Link](https://bedtools.readthedocs.io/en/latest/content/installation.html)] (>= v2.25.0)
+ Rscript [[Link](https://cran.r-project.org/mirrors.html)] (>= v3.2.0)
+ BWA [[Link](https://sourceforge.net/projects/bio-bwa/files)] (>=0.7.17)
+ Java (in order to run Picard and GATK)
+ Picard jar file [[Link](https://github.com/broadinstitute/picard)] (>= v2.23)
+ GATK jar file [[Link](https://github.com/broadinstitute/gatk/releases)] (>= v4.0.0.0, check compatibility with version of Java here: [[LINK](https://gatk.broadinstitute.org/hc/en-us/articles/360035532332-Java-version-issues)])

## Installation
The tool does not require compilation. It will download reference genome files from UCSC during the first run.

## Usage
This tool is made to be run by batch of samples ideally sequenced in the same batch or at least with same capture kit and sequencer. This is required for proper comparison between samples to detect deletions.
The sequencing data can be WES, WGS or targeted sequencing (if opsin locus, PRPF3 and PRPS1 are covered).

The main script 01_opsins-main.sh takes as input a tsv file with 2 or 3 columns: ID, forward FASTQ file and reverse FASTQ file (for paired-end). From these FASTQ files, Opsin-Miner will determine possible deletions, detect pathogenic haplotypes and coding variants.
It is called with bash:
```
bash 01_opsins-main.sh --script script_folder --out output_folder --input input_file.tsv --config config.txt
```
The approximate computation time per sample is few minutes per sample depending on the size of FASTQ files.

#### Required arguments
Option | Value  | Description 
--- | --- | ---
--script | STRING | Folder where the scripts are (e.g. XYZ/Opsin-Miner-main)
--out | STRING | Folder for outputs
--input | STRING | Tab-delimited text file with 2 or 3 columns: ID, forward FASTQ file and eventually reverse FASTQ file
--config | STRING | config text file containging path to executables

#### Optional argument
Option  | Description 
--- | ---
--malesonly | To be added if only males are analyzed (to avoid failing of sex identification).

## Outputs
The main output files which contain all required information are:
+ 0.events.tsv: list of detected events which can be rare or pathogenic variants, pathogenic haplotypes or deletions.
+ 0.phenotypes-all.tsv: inferred phenotype per sample based on detected event(s) as defined in “Determination of phenotype”.
+ 0.sex.tsv: inferred sex and Z-score of chrX coverage vs autosome (high = female, low = male)
+ 0.copy-number.tsv: number of copies of OPN1LW, OPN1MW and hybrid genes

The other output files are placed into multiple directories:
+ 00_raw-sequences: contains the selected sequences for each sample in FASTQ format
+ 01_bam: contains the aligned selected sequences for each sample in BAM format with indexes
+ 02_vcf: detected variants in VCF format for each sample
+ 03_variants: list of variants with information about reads for each sample and overall
+ 04_coverage: coverage for all regions for each sample and overall
+ 05_haplotypes: detected haplotypes (DNA and protein) for each sample and overall
+ 06_copy-number: dnumber of copies detected and plots with probability for each state
+ 07_plots: plots sumarizing information for each sample as well as Z-score for sex determination
+ 08_logs: logs of BWA, GATK and picard commands used

## Determination of phenotype

Phenotype is defined as “BCM” if there is:
+ Pathogenic haplotype in more of 75% of reads
+ Known pathogenic or LoF variant in more of 75% of reads
+ Deletion of LCR (Locus Control Region)

Phenotype is defined as “Color vision defiency or BCM“ if there is:
+ Deletions in LW and MW genes, a hybrid gene and a VUS in more of 75% of reads

Phenotype is defined as “Color vision defiency suggested” if there is:
+ Deletions in LW and MW genes and a hybrid gene without futher events
+ Deletions in LW and MW genes without futher events
+ Deletions in LW or MW genes

Phenotype is defined as “Inconclusive” if there is:
+ Pathogenic haplotype in less than 75% of reads
+ Known pathogenic or LoF variant in less than 75% of reads

Phenotype is defined as “Normal / Color vision deficency possible” if there is:
+ A rare variant (VUS) which is not know to be pathogenic in more than 10% of reads

Phenotype is defined as “Normal with VUS” if there is:
+ A rare variant (VUS) which is not know to be pathogenic in less than 10% of reads

List of pathogenic variants relative to OPN1LW (excluding LoF):
+ NP_064445.2:p.(Asn94Lys)
+ NP_064445.2:p.(Trp177Arg)
+ NP_064445.2:p.(Cys203Arg)
+ NP_064445.2:p.(Arg330Gln)
+ NP_064445.2:p.(Gly338Glu)


List of abreviations used:
+ BCM: Blue Cone Monochromacy
+ AF: Allele Frequency
+ VUS: Variant of Uncertain Significance
+ LoF: Loff-of-function

