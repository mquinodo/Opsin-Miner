# Opsin-Miner
Opsin-Miner is a bioinformatic tool that is used to analyze genotype at the OPN1LW/OPN1MW genes (red and green opsin) locus on the X chromosome. It identifies the number of gene copies, deletions, small variants, deletions, hybrid genes, and haplotypes and can be used with data from short-read targeted, exome, or genome sequencing.

This software was written by Mathieu Quinodoz in the group of Prof. Rivolta from the IOB in Basel, Switzerland in collaboration with Dr. Susanne Kohl and Prof. Bernd Wissinger from the University of Tübingen, Germany and mutliple other partners.

It was developed on Ubuntu 20.04.5 LTS (GNU/Linux 5.4.0-131-generic x86_64).

## Prerequisites (version recommended)
+ SAMtools [[Link](http://www.htslib.org/download/)] (>= v1.10)
+ BCFTools [[Link](https://samtools.github.io/bcftools/howtos/install.html)] (>= v1.9)
+ BEDTools [[Link](https://bedtools.readthedocs.io/en/latest/content/installation.html)] (>= v2.25.0)
+ Rscript [[Link](https://cran.r-project.org/mirrors.html)] (>= v3.2.0)
+ BWA [[Link](https://sourceforge.net/projects/bio-bwa/files)] (>=0.7.17)
+ Java (in order to run Picard and GATK)
+ Picard jar file [[Link](https://github.com/broadinstitute/picard)] (>= v2.23)
+ GATK jar file [[Link](https://github.com/broadinstitute/gatk/releases)] (>= v4.0.0.0, check compatibility with the version of Java here: [[LINK](https://gatk.broadinstitute.org/hc/en-us/articles/360035532332-Java-version-issues)])

## Installation
The tool does not require compilation. It will download reference files during the first run.

## Usage
This tool is made to be run by batch of samples ideally sequenced in the same batch or at least with the same capture kit and sequencer. For optimal results, we recommend analyzing at least 10 individuals, either all males (--malesonly option) or a balanced cohort including at least three males and three females.

The main script 01_opsins-main.sh takes as input a tsv file with 2 or 3 columns: ID, BAM/CRAM file or forward FASTQ file, and reverse FASTQ file (for paired-end). 
It is called with bash:
```
bash 01_opsins-main.sh --script script_folder --out output_folder --input input_file.tsv --config config.txt
```
The approximate computation time per sample is few seconds per sample for BAM/CRAM files and few minutes for FASTQ files.

#### Required arguments
Option | Value  | Description 
--- | --- | ---
--script | STRING | Folder where the scripts are (e.g. /home/user/Opsin-Miner-main)
--out | STRING | Folder for outputs
--input | STRING | Tab-delimited text file with 2 or 3 columns: ID, BAM/CRAM file or forward FASTQ file and eventually reverse FASTQ file
--config | STRING | Configuration text file containing path to executables
--fastq | STRING | To be used for the analysis FASTQ files (use either --fastq or --bam-cram)
--bam-cram | STRING | To be used for the analysis BAM or CRAM files (use either --fastq or --bam-cram)

#### Optional important arguments
Option | Value  | Description 
--- | --- | ---
--malesonly | NA | To be added if only males are analyzed (to avoid failure of sex identification).
--hg38 | NA | To be used if BAM/CRAM files are used and mapped to hg38 (and not hg19)
--WGS | NA | To be used with data from genome sequencing (as opposed to exome or targeted sequencing)
--noCHR | NA | To be used if BAM/CRAM files have a chromosome notation without "chr" ("1" instead of "chr1")

#### Optional arguments
Option | Value  | Description 
--- | --- | ---
--smallLCR | NA | To be used with exome or targeted sequencing data if the LCR is captured, in order to detect smaller deletions.
--sexfile | STRING | Text files with sex provided if known with certainty for all samples (tab-delimited with ID and Male/Female).
--normauto | STRING | If chromosome 1 cannot be used for normalization, a gene name on an autosome to be used (not recommended).
--normX | STRING | If chromosome X cannot be used for normalization, a gene name on chromosome X to be used (not recommended).

## Outputs
The main output files which contain all required information are:
+ 0.events.tsv: list of detected events which can be rare or pathogenic variants, pathogenic haplotypes or deletions.
+ 0.phenotypes.tsv: inferred phenotype per sample based on detected event(s) as defined in “Determination of phenotype”.

The other output files are placed into multiple directories:
+ 00_raw-sequences: contains the selected sequences for each sample in FASTQ format
+ 01_bam: contains the realigned selected sequences for each sample in BAM format with indexes (hg19). The OPN1LWMW BAM is aligned on the modified hg19 reference with masked OPN1MW2.
+ 02_vcf: detected variants in VCF format for each sample
+ 03_variants: list of variants with information about reads for each sample and overall
+ 04_coverage: coverage for all regions for each sample and overall as well as Z-score for sex determination
+ 05_haplotypes: detected haplotypes (DNA and protein) for each sample and overall
+ 06_copy-number: number of copies detected and plots with probability for each state
+ 07_plots: plots summarizing information for each sample
+ 08_logs: logs of bcftools, BWA, GATK and picard commands used

List of abreviations used:
+ BCM: Blue Cone Monochromacy
+ BED: Bornholm Eye Disease
+ CVD: Color Vision Deficiency
+ AF: Allele Frequency
+ VUS: Variant of Uncertain Significance
+ LoF: Loss-of-function

