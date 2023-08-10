#!/bin/bash

usage() { echo "## ERROR: Usage: $0 [--script <string>] [--out <string>] [--input <string>] [--config <string>] [--malesonly]. Exit." 1>&2; exit 1; }

while getopts ":-:" o; do
    case "${o}" in
    	-)  
        	case $OPTARG in
                script)
                    script="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                out)
                    out="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                input)
                    input="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                config)
                    config="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
               	malesonly)
                    malesonly="Yes"
                    ;;
				*)
            		usage
            	;;
         esac ;;        
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${script}" ]; then
    usage
fi
if [ -z "${out}" ]; then
    usage
fi
if [ -z "${input}" ]; then
    usage
fi
if [ -z "${config}" ]; then
    usage
fi

if [ ! -f "$config" ]; then
  echo "Config file does not exist!"
  exit 1
fi

if [ -z "${malesonly}" ]; then
    echo " -> Analysis will be made assume all samples are males."
    malesonly="No"
fi

samtools=$(grep "^samtools" "$config" | cut -f2)
bcftools=$(grep "^bcftools" "$config" | cut -f2)
bedtools=$(grep "^bedtools" "$config" | cut -f2)
Rscript=$(grep "^Rscript" "$config" | cut -f2)
bwa=$(grep "^bwa" "$config" | cut -f2)
java=$(grep "^java" "$config" | cut -f2)
picard=$(grep "^picard" "$config" | cut -f2)
gatk=$(grep "^gatk" "$config" | cut -f2)

nosamtools() { echo "## ERROR: samtools lower than v1.10 -> Please Update! Exit." 1>&2; exit 1; }
nobcftools() { echo "## ERROR: bcftools lower than v1.9 -> Please Update! Exit." 1>&2; exit 1; }
nobedtools() { echo "## ERROR: bedtools lower than v2.24.0 -> Please Update! Exit." 1>&2; exit 1; }
noRversion() { echo "## ERROR: R lower than v3.2.0 -> Please Update! Exit." 1>&2; exit 1; }
nobwa() { echo "## ERROR: BWA not installed -> Please Install! Exit." 1>&2; exit 1; }
nojava() { echo "## ERROR: Java not installed -> Please Install! Exit." 1>&2; exit 1; }
nopicard() { echo "## ERROR: picard jar file not found -> Please download it! Exit." 1>&2; exit 1; }
nogatk() { echo "## ERROR: gatk jar file not found -> Please download it Exit." 1>&2; exit 1; }

currentver="$($samtools --version | head -n1 | cut -d" " -f2)"
requiredver="1.10"
if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then 
    echo "# samtools higher or equal to v1.10"
else
    nosamtools
fi

currentver="$($bcftools -v | head -n1 | cut -d" " -f2)"
requiredver="1.9"
if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then 
    echo "# bcftools higher or equal to v1.9"
else
    nobcftools
fi

currentver="$($bedtools --version | cut -d" " -f2)"
requiredver="v2.24.0"
if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then 
    echo "# bedtools higher or equal to v2.24.0"
else
    nobedtools
fi

currentver="$($Rscript --version | grep "version" | cut -d" " -f4)"
requiredver="3.2.0"
if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then 
    echo "# R higher or equal to v3.2.0"
else
    noRversion
fi

if ! command -v $bwa &> /dev/null
then
    nobwa
else
	echo "# BWA installed"
fi

if ! command -v $java &> /dev/null
then
    nojava
else
	echo "# java installed"
fi

if [ ! -f "$picard" ]; then
    nopicard
else
	echo "# picard jar file found"
fi

if [ ! -f "$gatk" ]; then
    nogatk
else
	echo "# gatk jar file found"
fi


mkdir -p $out
mkdir -p $out/00_raw-sequences
mkdir -p $out/01_bam
mkdir -p $out/02_vcf
mkdir -p $out/03_variants
mkdir -p $out/04_coverage
mkdir -p $out/05_haplotypes
mkdir -p $out/06_copy-number
mkdir -p $out/07_plots
mkdir -p $out/08_logs
mkdir -p $out/temp

if [ ! -f $script/data/hg19.p13.plusMT.no_alt_analysis_set.fa ]
then
	echo "Step 0.1: Downloading reference genome (will be done only once)"
	wget -q -P $script/data https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz
	gunzip $script/data/hg19.p13.plusMT.no_alt_analysis_set.fa.gz
fi

if [ ! -f $script/data/hg19.p13.plusMT.no_alt_analysis_set.fa.amb ]
then
	echo "Step 0.2: Downloading BWA index (will be done only once)"
	wget -q -P $script/data https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.bwa_index.tar.gz
	tar -xf $script/data/hg19.p13.plusMT.no_alt_analysis_set.bwa_index.tar.gz -C $script/data/
	mv $script/data/hg19.p13.plusMT.no_alt_analysis_set/* $script/data/
	rm -f $script/data/hg19.p13.plusMT.no_alt_analysis_set.bwa_index.tar.gz
	rm -rf $script/data/hg19.p13.plusMT.no_alt_analysis_set
fi

ref1=$script/data/hg19.p13.plusMT.no_alt_analysis_set.fa
ref2=$script/data/chrX.masked-MW2.fa


if [ ! -f $ref2 ]
then
	echo "Step 0.3: Creating chrX reference with masked OPN1MW2 and polymorphism (will be done only once)"
	# extract chrX
	$samtools faidx $ref1 chrX > $script/data/chrX.fa

	# mask region with OPN1MW and OPN1MW2
	$bedtools maskfasta -fi $script/data/chrX.fa -bed $script/data/OPN1MW2.bed -fo $script/data/chrX.temp.fa

	# remplace difference in exon 5
	sed 's/CCATACGCCTTCTTCGCATGCTTTGCTGCTGCCAACCCTGGCTACCCCTTCCACCCTTTG/CCNTACGCCTTCTTCGCATGCTTTGCTGCTGCCAACCCTGGCTACCCCTTCCACCCTTTG/g' $script/data/chrX.temp.fa > $script/data/chrX.masked-MW2.fa

	# remove temp files
	rm $script/data/chrX.temp* $script/data/chrX.fa
	# indexes for the chrX modified reference
	$samtools faidx $script/data/chrX.masked-MW2.fa
	$bwa index $script/data/chrX.masked-MW2.fa 2>> $out/08_logs/BWA-index-genome.txt
	$java -jar $picard CreateSequenceDictionary -R $script/data/chrX.masked-MW2.fa -O $script/data/chrX.masked-MW2.dict >> $out/08_logs/picard-dictionary.txt 2>&1
fi


IFS=$'\r\n' GLOBIGNORE='*' command eval 'patients=($(cut -f1 $input))'
IFS=$'\r\n' GLOBIGNORE='*' command eval 'fastq1=($(cut -f2 $input))'
IFS=$'\r\n' GLOBIGNORE='*' command eval 'fastq2=($(cut -f3 $input))'
n=$((${#patients[@]}-1))

echo "Step 1: Greping sequences in FASTQ files"
for i in `seq 0 $n`
do
	pat=${patients[$i]}

	echo "   Processing sample: $pat"
	# greping all sequences
	if [ ! -f $out/00_raw-sequences/$pat.sel.fastq.gz ]
	then
		zcat ${fastq1[$i]} | grep -B1 -A2 -Ff $script/data/region.grep.fa > $out/00_raw-sequences/${pat}_1.tsv &
		pid1=$!
		if [ ! -z "${fastq2[$i]}" ]
		then
			zcat ${fastq2[$i]} | grep -B1 -A2 -Ff $script/data/region.grep.fa > $out/00_raw-sequences/${pat}_2.tsv &
			pid2=$!
			wait $pid2
		fi
		wait $pid1
		
		if [ ! -z "${fastq2[$i]}" ]
		then
			cat $out/00_raw-sequences/${pat}_1.tsv $out/00_raw-sequences/${pat}_2.tsv > $out/00_raw-sequences/$pat.sel.fastq
		else
			cat $out/00_raw-sequences/${pat}_1.tsv > $out/00_raw-sequences/$pat.sel.fastq
		fi
		
		gzip $out/00_raw-sequences/$pat.sel.fastq
		rm -f $out/00_raw-sequences/${pat}_1.tsv $out/00_raw-sequences/${pat}_2.tsv

	fi

done

echo "Step 2: Determining coverage and haplotypes"
for pat in "${patients[@]}"
do
	echo "   Processing sample: $pat"

	file=$out/01_bam/$pat.bam
	
	if [ ! -f $file ]
	then
		# recreate BAM file: align selected reads and then convert SAM to BAM
		RG=$(echo "@RG\tID:$pat\tSM:$pat\tPL:Illumina")
		$bwa mem -M -C -t 32 -R $RG $ref1 $out/00_raw-sequences/$pat.sel.fastq.gz > $out/temp/$pat.sam 2>> $out/08_logs/BWA-genome.txt
		$java -jar $picard SortSam -VALIDATION_STRINGENCY SILENT -SORT_ORDER coordinate -I $out/temp/$pat.sam -O $out/01_bam/$pat.bam -TMP_DIR $out/temp >> $out/08_logs/picard-genome.txt 2>&1
		$samtools index $out/01_bam/$pat.bam
		rm $out/temp/$pat.sam
	fi

	if [ ! -f $out/05_haplotypes/$pat.haplotypes.tsv ]
	then
		# extract sequences mapping to OPN1 region
		$samtools view $out/01_bam/$pat.bam "chrX:153390252-153510581" | cut -f10 > $out/temp/$pat.seq.all.tsv

		# exon3 OPN1LW sense based on two regions with changes in haplotypes
		grep -P "GTGAGATTTGATGCCAAGCTGGCCATC|TGCAAGCCCTTTGGCAA" $out/temp/$pat.seq.all.tsv > $out/temp/$pat.seq.hap-fwd.tsv
		# exon3 OPN1LW antisense
		grep -P "GATGGCCAGCTTGGCATCAAATCTCAC|TTGCCAAAGGGCTTGCA" $out/temp/$pat.seq.all.tsv | tr ACGTacgt TGCAtgca | rev > $out/temp/$pat.seq.hap-rev.tsv
		# merge sense and antisense
		cat $out/temp/$pat.seq.hap-fwd.tsv $out/temp/$pat.seq.hap-rev.tsv | sort | uniq > $out/temp/$pat.seq.hap-all.tsv
		# grep sequence that are long enough on both sides, first left and then right
		grep -P "GTGGCTGGTGGT|ATGGATGGTGGT|GTGGATGGTGGT|ATGGCTGGTGGT" $out/temp/$pat.seq.hap-all.tsv | grep -P "CTTCTCCTGGATCTGGT|CTTCTCCTGGGTCTGGT|CTTCTCCTGGGTCTGGG|CTTCTCCTGGATCTGGG" > $out/temp/$pat.sequences.tsv
		# R script for sequence alignement and haplotype calling
		$Rscript $script/sequence-analysis.R $out/temp/$pat.sequences.tsv $out/05_haplotypes/$pat.haplotypes.tsv $out/05_haplotypes/$pat.haplotypesDNA.tsv
		# remove temporary files
		rm $out/temp/$pat.seq.hap*
	fi

	if [ ! -f $out/04_coverage/$pat.coverage.tsv ]
	then
		# look at coverage
		cat $out/temp/$pat.seq.all.tsv | tr ACGTacgt TGCAtgca | rev > $out/temp/$pat.seq.all-rev.tsv
		cat $out/temp/$pat.seq.all.tsv $out/temp/$pat.seq.all-rev.tsv > $out/temp/$pat.seq.all-both.tsv

		grep CACAGGCCAGTATAAAGCGCCGTGAC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon1-part1" "\t" $1}' > $out/04_coverage/$pat.coverage.tsv
		grep CACTGGCCGGTATAAAGCACCGTGAC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon1-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CTCAGGTGATGCGCCAGGGCCGGCT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon1-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CTCAGGTGACGCACCAGGGCCGGCT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon1-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep GATGATCTTTGTGGTCACTGCATCCGTC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon2-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep GATGATCTTTGTGGTCATTGCATCCGTC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon2-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep GTGAACCAGGTCTCTGGCTACTTCGTGCTGGG $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon2-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep GTGAACCAGGTCTATGGCTACTTCGTGCTGGG $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon2-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep GTGAGATTTGATGCCAAGCTGGCCATC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1-exon3-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep TGTGTGGACAGCCCCGCCCATCTTTGGTTGGAGCAG $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1-exon3-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep CCTGCTGCATCATCCCACTCGCT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon4-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CCTGCTGCATCACCCCACTCAGC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon4-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep GCTATCATCATGCTCTGCTACCT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon4-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep AGCATCATCGTGCTCTGCTACCT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon4-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep TACTGCGTCTGCTGGGGACCC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon5-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep TTCTGCTTCTGCTGGGGACCA $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MW-exon5-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep TTCTGCTTCTGCTGGGGACCC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MW2-exon5-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CTACACCTTCTTCGCATGCTT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon5-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep ATACGCCTTCTTCGCATGCTT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MW-exon5-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CTACGCCTTCTTCGCATGCTT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MW2-exon5-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep GAAACTGCATCTTGCAGCTTT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1-exon6-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep AGGTCTCATCTGTGTCCTCGG $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1-exon6-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		$samtools view $out/01_bam/$pat.bam "chrX:106869654-106896256" | wc -l | awk -F"\t" '{print "Control-chrX-PRPS1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		$samtools view $out/01_bam/$pat.bam "chr1:150291917-150327704" | wc -l | awk -F"\t" '{print "Control-autosome-PRPF3" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		$samtools view $out/01_bam/$pat.bam "chrX:153404000-153408000" | wc -l | awk -F"\t" '{print "LCR" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		$samtools view $out/01_bam/$pat.bam "chrX:143653980-143658741" | wc -l | awk -F"\t" '{print "LCR-control" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

	fi
	if [ ! -f $out/04_coverage/$pat.coverage.mut.tsv ]
	then
		$samtools view $out/01_bam/$pat.bam "chrX:153390252-153510581" | cut -f10 > $out/temp/$pat.seq.all.tsv
		cat $out/temp/$pat.seq.all.tsv | tr ACGTacgt TGCAtgca | rev > $out/temp/$pat.seq.all-rev.tsv
		cat $out/temp/$pat.seq.all.tsv $out/temp/$pat.seq.all-rev.tsv > $out/temp/$pat.seq.all-both.tsv
		# mutations
		grep CGCGGCCCAGACGTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "C203R-LW" "\t" $1}' > $out/04_coverage/$pat.coverage.mut.tsv
		grep CGCGGCCCAGACGTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "C203R-MW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		grep TGCGGCCCAGACGTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "C203-LW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		grep TGCGGCCCAGACGTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "C203-MW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		grep ATGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "V207M-LW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		grep ATGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "V207M-MW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		grep GTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAT $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "V207-LW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		grep GTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "V207-MW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		#grep ATGGTCATCTGCTGCATCATCCCACTC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "T226I-LW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		#grep ATGGTCATCTGCTGCATCACCCCACTC $out/temp/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "T226I-MW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
	fi

	rm -f $out/temp/$pat.*
	
done

echo "Step 3: Analysis of rare variants"
for pat in "${patients[@]}"
do
	echo "   Processing sample: $pat"
	if [ ! -f $out/03_variants/$pat.variants.tsv ]
	then
		# take reads from BAM (OPN1 region) and create FASTQs
		$samtools view -b -h $out/01_bam/$pat.bam "chrX:153390252-153510581" > $out/temp/$pat.region.bam
		$samtools index $out/temp/$pat.region.bam
		$samtools sort $out/temp/$pat.region.bam > $out/temp/$pat.region.sorted.bam
		$samtools index $out/temp/$pat.region.sorted.bam
		$bedtools bamtofastq -i $out/temp/$pat.region.sorted.bam -fq $out/temp/$pat.region.fastq

		# alignement, sam to bam and indexing
		RG=$(echo "@RG\tID:$pat\tSM:$pat\tPL:Illumina")
		$bwa mem -M -C -t 32 -R $RG $ref2 $out/temp/$pat.region.fastq > $out/temp/$pat.region.sam 2>> $out/08_logs/BWA-region.txt
		$java -jar $picard SortSam -VALIDATION_STRINGENCY SILENT -SORT_ORDER coordinate -I $out/temp/$pat.region.sam -O $out/temp/$pat.region2.bam -TMP_DIR $out/temp 2>> $out/08_logs/picard-region.txt
		$samtools index $out/temp/$pat.region2.bam

		# variant calling
		$java -jar $gatk HaplotypeCaller \
		-I $out/temp/$pat.region2.bam \
		-R $ref2 \
		-O $out/temp/$pat.vcf \
		-G StandardAnnotation \
		-G StandardHCAnnotation \
		-L $script/data/hg19_ncbiRefSeq-exons.LW-MW.bed 2>> $out/08_logs/HaplotypeCaller.txt

		# VCF processing
		$bcftools norm -m -both --fasta-ref $ref2 -o $out/02_vcf/$pat.norm.vcf $out/temp/$pat.vcf 2>> $out/08_logs/bcftools.txt
		grep -P "chrX\t" $out/02_vcf/$pat.norm.vcf | awk -F"\t" -v pat="$pat" '{split($10,a,":"); split(a[2],b,","); print pat "\t" $1 "\t" $2 "\t" $4 "\t" $5 "\t" b[1] "\t" b[2] "\t" $1 "-" $2 "-" $4 "-" $5}' > $out/03_variants/$pat.variants.tsv
		rm -rf $out/temp/$pat.*

	fi

done

echo "Step 4: Batch processing"

# merge variants in one file
echo -e "ID\tchr\tpos\tref\talt\treads-ref\treads-alt\tchange" > $out/03_variants/0.variants-all.tsv
cat $out/03_variants/*.variants.tsv >> $out/03_variants/0.variants-all.tsv

echo -e "ID\tchr\tpos\tref\talt\treads-ref\treads-alt\tchange" > $out/03_variants/0.variants-rare.tsv
cat $out/03_variants/*.variants.tsv | grep -v -Ff $script/data/gnomAD-frequent.tsv >> $out/03_variants/0.variants-rare.tsv

cat $out/04_coverage/*.coverage.mut.tsv > $out/04_coverage/0.coverage.mut.ALL.tsv

# extract pathogenic and rare haplotypes
grep -P "A|I|V|S" $out/05_haplotypes/*.haplotypes.tsv | awk -F"\t" '{split($1,a,":"); n=split(a[1],b,"/"); split(b[n],c,"."); print c[1] "\t" a[2] "\t" $2}' > $out/05_haplotypes/0.PROT.all.tsv
grep -P "LIAVA|LVAVA|MIAVA|LIAVS|LIVVA" $out/05_haplotypes/0.PROT.all.tsv > $out/05_haplotypes/0.PROT.pathogenic.tsv
echo -e "NA\tNA\t0" >> $out/05_haplotypes/0.PROT.pathogenic.tsv
cat $out/05_haplotypes/*.haplotypes.tsv | cut -f1 | sort | uniq -c | gawk '$1<=4{print $2}' > $out/temp/0.PROT.temp.tsv
grep -Ff $out/temp/0.PROT.temp.tsv $out/05_haplotypes/0.PROT.all.tsv | awk -F"\t" '{if($3>4) print $0}' > $out/05_haplotypes/0.PROT.rare.tsv

# extract rare DNA haplotypes
grep -P "T|C|G" $out/05_haplotypes/*.haplotypesDNA.tsv | awk -F"\t" '{split($1,a,":"); n=split(a[1],b,"/"); split(b[n],c,"."); print c[1] "\t" a[2] "\t" $2}' > $out/05_haplotypes/0.DNA.all.tsv
cat $out/05_haplotypes/*.haplotypesDNA.tsv | cut -f1 | sort | uniq -c | gawk '$1<=4{print $2}' > $out/temp/0.DNA.temp.tsv
grep -Ff $out/temp/0.DNA.temp.tsv $out/05_haplotypes/0.DNA.all.tsv | awk -F"\t" '{if($3>4) print $0}' > $out/05_haplotypes/0.DNA.rare.tsv

# merge coverage files
printf "Name" > $out/temp/header.tsv
cut -f1 $out/04_coverage/${patients[0]}.coverage.tsv > $out/temp/temp0.tsv
cat $out/temp/temp0.tsv > $out/temp/data.tsv
for pat in "${patients[@]}"
do
	cut -f2 $out/04_coverage/$pat.coverage.tsv > $out/temp/temp.tsv
	cat $out/temp/data.tsv > $out/temp/temp2.tsv
	paste -d"\t" $out/temp/temp2.tsv $out/temp/temp.tsv > $out/temp/data.tsv
	printf "\t$pat" >> $out/temp/header.tsv
done
echo "" >> $out/temp/header.tsv
sort -k1,1 $out/temp/data.tsv > $out/temp/data2.tsv
cat $out/temp/header.tsv $out/temp/data2.tsv > $out/04_coverage/0.coverage-ALL.tsv

echo "Step 5: Creating plots and events files"
$Rscript $script/analysis_plots.R $out/04_coverage/0.coverage-ALL.tsv $out/05_haplotypes/0.PROT.all.tsv $out/05_haplotypes/0.PROT.pathogenic.tsv $out/05_haplotypes/0.DNA.all.tsv $out/03_variants/0.variants-rare.tsv $script/data/HGVS-gnomAD.RData $out $out/04_coverage/0.coverage.mut.ALL.tsv $malesonly

rm -rf $out/temp

# rs949430 rs713 rs731614 rs148267962 rs5986963 rs5986964 rs149897670 rs145009674 rs949431
# R151R L153M V155V N161N V171M V171V A174V I178V S180A
