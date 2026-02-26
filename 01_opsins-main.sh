#!/bin/bash

usage() { echo "## ERROR: Usage: $0 [--script <string>] [--out <string>] [--input <string>] [--config <string>] [--hg38] [--normauto <string>] [--normX <string>] [--malesonly] [--fastq] [--bam-cram] [--smallLCR] [--WGS] [--noCHR] [--sexfile]. Exit." 1>&2; exit 1; }

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
                normauto)
                    normauto="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                normX)
                    normX="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                sexfile)
                    sexfile="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                hg38)
                    hg38="Yes"
                    ;;
               	malesonly)
                    malesonly="Yes"
                    ;;
                bam-cram)
                    bam="Yes"
                    ;;
                fastq)
                    fastq="Yes"
                    ;;
                smallLCR)
                    smallLCR="Yes"
                    ;;
                WGS)
                    WGS="Yes"
                    ;;
				noCHR)
                    noCHR="Yes"
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

if [ ! -z "${hg38}" ]; then
    echo " -> BAM / CRAM files aligned on hg38."
fi

if [ ! -z "${sexfile}" ]; then
    echo " -> Sex will be taken from provided file."
fi
if [ -z "${sexfile}" ]; then
    echo " -> Sex will be infered from sequencing data."
    sexfile="No"
fi

if [ ! -z "${normauto}" ] && [ ! -z "${fastq}" ]; then
    echo " -> Normalization on autosomes will be done using custom gene."
fi
if [ ! -z "${normauto}" ] && [ -z "${fastq}" ]; then
    echo " -> Normalization on autosomes will be done using custom gene."
fi
if [ -z "${normauto}" ] && [ ! -z "${fastq}" ]; then
    echo " -> Normalization on autosomes will be done using PRPF31 gene."
    normauto="PRPF31"
fi
if [ -z "${normauto}" ] && [ -z "${fastq}" ]; then
    echo " -> Normalization on autosomes will be done using full chromosome 1."
    normauto="chr1"
fi

if [ ! -z "${normX}" ] && [ ! -z "${fastq}" ]; then
    echo " -> Normalization on chromosome X will be done using custom gene."
fi
if [ ! -z "${normX}" ] && [ -z "${fastq}" ]; then
    echo " -> Normalization on chromosome X will be done using custom gene."
fi
if [ -z "${normX}" ] && [ ! -z "${fastq}" ]; then
    echo " -> Normalization on chromosome X will be done using CACNA1F gene."
    normX="CACNA1F"
fi
if [ -z "${normX}" ] && [ -z "${fastq}" ]; then
    echo " -> Normalization on chromosome X will be done using full chromosome X."
    normX="chrX"
fi

if [ ! -z "${malesonly}" ]; then
    echo " -> Analysis will be made assuming all samples are males."
fi
if [ -z "${malesonly}" ]; then
    malesonly="No"
fi

if [ ! -z "${bam}" ]; then
    echo " -> Analysis will be made with BAM or CRAM files."
fi

if [ -z "${fastq}" ] && [ -z "${bam}" ]; then
    echo "## ERROR: Please indicate if using FASTQ, BAM or CRAM files. Exit." 1>&2; exit 1; 
fi

if [ ! -z "${fastq}" ]; then
    echo " -> Analysis will be made with FASTQ files."
    bam="No"
fi

if [ ! -z "${smallLCR}" ]; then
    region3="chrX:153405983-153406393"
fi
if [ -z "${smallLCR}" ]; then
    region3="chrX:153404000-153408000"
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

if [ ! -f $script/data/bbmap/bbmerge.sh ]
then
	echo "Step 0.3: Downloading BBMap (will be done only once)"
	wget -q -P $script/data https://storage.googleapis.com/dominofiles/BBMap_39.06.tar.gz
	tar -xzf $script/data/BBMap_39.06.tar.gz -C $script/data
	rm -rf $script/data/BBMap_39.06.tar.gz
fi

ref1=$script/data/hg19.p13.plusMT.no_alt_analysis_set.fa
ref2=$script/data/chrX.masked-MW2.fa

if [ ! -f $ref2 ]
then
	echo "Step 0.4: Creating chrX reference with masked OPN1MW2 and polymorphism (will be done only once)"
	# extract chrX
	$samtools faidx $ref1 chrX > $script/data/chrX.fa

	# mask region with OPN1MW and OPN1MW2
	$bedtools maskfasta -fi $script/data/chrX.fa -bed $script/data/OPN1MW2.bed -fo $script/data/chrX.temp.fa

	# replace difference in exon 5 by N to match both OPN1MW and OPN1MW2
	sed 's/CCATACGCCTTCTTCGCATGCTTTGCTGCTGCCAACCCTGGCTACCCCTTCCACCCTTTG/CCNTACGCCTTCTTCGCATGCTTTGCTGCTGCCAACCCTGGCTACCCCTTCCACCCTTTG/g' $script/data/chrX.temp.fa > $script/data/chrX.masked-MW2.fa

	# remove temp files
	rm $script/data/chrX.temp* $script/data/chrX.fa
	# indexes for the chrX modified reference
	$samtools faidx $script/data/chrX.masked-MW2.fa
	$bwa index $script/data/chrX.masked-MW2.fa 2>> $out/08_logs/BWA-index-genome.txt
	$java -jar $picard CreateSequenceDictionary -R $script/data/chrX.masked-MW2.fa -O $script/data/chrX.masked-MW2.dict >> $out/08_logs/picard-dictionary.txt 2>&1
fi

echo "Step 1: Creating fasta file for sequences to be searched in reads"
# preparing file with sequences from opsin cluster, one autosomal control gene and one chrX control gene
if [ ! -f "$script/data/hg19_ncbiRefSeq.tsv" ]; then
    echo "Unzipping of the gene file."
    gunzip -q $script/data/hg19_ncbiRefSeq.tsv.gz
fi


if [ ! -z "${hg38}" ] && [ ! -f "$out/00_raw-sequences/0.region.hg38.${normauto}.${normX}.bed" ]; then
	if [ ! -f "$script/data/hg38_ncbiRefSeq.tsv" ]; then
    	echo "Unzipping of the gene file."
    	gunzip -q $script/data/hg38_ncbiRefSeq.tsv.gz
	fi
    file=$out/temp/hg38_ncbiRefSeq
    initial=$script/data/hg38_ncbiRefSeq.tsv
	grep "NM_" $initial | awk -F"\t" '{n=split($10,start,","); split($11,end,","); for(i=1; i<n; i=i+1) {ex=i; if($4=="-"){ex=(n-i)}; if($3 !~ /_/) {print $3 "\t" start[i]-00 "\t" end[i]+00 "\t" $13 "_" $2 "_exon" ex}  } }' | awk -F"\t" '{if($3>$2) print $0}' > $file.exons.bed

	# large regions for WGS
	if [ ! -z "${WGS}" ]; then
		echo -e "chrX\t154137280\t154160853\tLCR-OPN1LW" > $file.exons.sel.bed # OPN1LW+LCR
		echo -e "chrX\t154178338\t154197722\tOPN1MW" >> $file.exons.sel.bed # OPN1MW
		echo -e "chrX\t154215399\t154234788\tOPN1MW2" >> $file.exons.sel.bed # OPN1MW2
	fi

	# small regions for WES
	if [ -z "${WGS}" ]; then
		# NM_020061.6 OPN1LW, NM_000513.2 OPN1MW, NM_001048181.3 OPN1MW2
		grep "NM_020061.6" $file.exons.bed > $file.exons.sel.bed
		grep "NM_000513.2" $file.exons.bed >> $file.exons.sel.bed
		grep "NM_001048181.3" $file.exons.bed >> $file.exons.sel.bed
		echo -e "chrX\t154138526\t154142526\tLCR" >> $file.exons.sel.bed # LCR
	fi

	if [ "${normauto}" != "chr1" ]; then
		grep "[^\t]${normauto}_" $file.exons.bed >> $file.exons.sel.bed
	fi
	if [ "${normX}" != "chrX" ]; then
		grep "[^\t]${normX}_" $file.exons.bed >> $file.exons.sel.bed
	fi

	# padding of 1000bp
	awk -F"\t" '{print $1 "\t" $2-1000 "\t" $3+1000 "\tNA"}' $file.exons.sel.bed > $file.exons.sel.pro2.bed
	cp $file.exons.sel.pro2.bed $out/00_raw-sequences/0.region.hg38.${normauto}.${normX}.bed
	rm -f $file.exons.*
fi

if [ ! -f "$out/00_raw-sequences/0.region.${normauto}.${normX}.bed" ]; then
	file=$out/temp/hg19_ncbiRefSeq
	initial=$script/data/hg19_ncbiRefSeq.tsv

	grep "NM_" $initial | awk -F"\t" '{n=split($10,start,","); split($11,end,","); for(i=1; i<n; i=i+1) {ex=i; if($4=="-"){ex=(n-i)}; if($3 !~ /_/) {print $3 "\t" start[i]-00 "\t" end[i]+00 "\t" $13 "_" $2 "_exon" ex}  } }' | awk -F"\t" '{if($3>$2) print $0}' > $file.exons.bed
	# large regions for WGS
	if [ ! -z "${WGS}" ]; then
		echo -e "chrX\t153402755\t153426328\tLCR-OPN1LW" > $file.exons.sel.bed # OPN1LW+LCR
		echo -e "chrX\t153443826\t153463213\tOPN1MW" >> $file.exons.sel.bed # OPN1MW
		echo -e "chrX\t153480869\t153500257\tOPN1MW2" >> $file.exons.sel.bed # OPN1MW2
	fi

	# small regions for WES
	if [ -z "${WGS}" ]; then
		# NM_020061.6 OPN1LW, NM_000513.2 OPN1MW, NM_001048181.3 OPN1MW2
		grep "NM_020061.6" $file.exons.bed > $file.exons.sel.bed
		grep "NM_000513.2" $file.exons.bed >> $file.exons.sel.bed
		grep "NM_001048181.3" $file.exons.bed >> $file.exons.sel.bed
		echo -e "chrX\t153404000\t153408000\tLCR" >> $file.exons.sel.bed # LCR
	fi

	# take exons from normalisation genes if chr1 and chrX (defaults) are not used
	if [ "${normauto}" != "chr1" ]; then
		grep "[^\t]${normauto}_" $file.exons.bed >> $file.exons.sel.bed
	fi
	if [ "${normX}" != "chrX" ]; then
		grep "[^\t]${normX}_" $file.exons.bed >> $file.exons.sel.bed
	fi

	# padding of 100bp
	awk -F"\t" '{print $1 ":" $2-100 "-" $3+100}' $file.exons.sel.bed > $file.exons.sel.pro.bed
	awk -F"\t" '{print $1 "\t" $2-1000 "\t" $3+1000 "\tNA"}' $file.exons.sel.bed > $file.exons.sel.pro2.bed

	rm -f $out/temp/region.grep.temp1.fa
	touch $out/temp/region.grep.temp1.fa
	while read p; do
	  samtools faidx $ref1 $p | grep -v ">" >> $out/temp/region.grep.temp1.fa
	done <$file.exons.sel.pro.bed

	# split in 30 bp regions
	fold -w 30 $out/temp/region.grep.temp1.fa | awk -F"\t" '{n=split($1,a,""); if(n==30) print $1}' > $out/temp/region.grep.temp2.fa

	# add exon 3 haplotypes
	echo "TCTCTGGCCATCATTTCCTGGGAGAG" >> $out/temp/region.grep.temp2.fa
	echo "TGGTGGTGTGCAAGCCCTTTGGCAA" >> $out/temp/region.grep.temp2.fa
	echo "TGGTGGTCTGCAAGCCCTTTGGCAA" >> $out/temp/region.grep.temp2.fa
	echo "GATTTGATGCCAAGCTGGCCATC" >> $out/temp/region.grep.temp2.fa
	echo "CTGCTGTGTGGACAGCCCCGCCCA" >> $out/temp/region.grep.temp2.fa

	cat $out/temp/region.grep.temp2.fa | tr ACGTacgt TGCAtgca | rev > $out/temp/region.grep.temp3.fa
	cat $out/temp/region.grep.temp2.fa $out/temp/region.grep.temp3.fa | sort | uniq > $out/00_raw-sequences/0.region.grep.${normauto}.${normX}.fa
	cp $file.exons.sel.pro2.bed $out/00_raw-sequences/0.region.${normauto}.${normX}.bed

	if [ "${normauto}" != "chr1" ]; then
		grep $normauto $file.exons.bed > $out/00_raw-sequences/0.region.exons.bed
	fi
	if [ "${normX}" != "chrX" ]; then
		grep $normX $file.exons.bed >> $out/00_raw-sequences/0.region.exons.bed
	fi
fi

# take regions from normalisation genes if chr1 and chrX (defaults) are not used
file=$out/temp/hg19_ncbiRefSeq
initial=$script/data/hg19_ncbiRefSeq.tsv
grep "NM_" $initial | awk -F"\t" '{n=split($10,start,","); split($11,end,","); for(i=1; i<n; i=i+1) {ex=i; if($4=="-"){ex=(n-i)}; if($3 !~ /_/) {print $3 "\t" start[i]-00 "\t" end[i]+00 "\t" $13 "_" $2 "_exon" ex}  } }' | awk -F"\t" '{if($3>$2) print $0}' > $file.exons.bed
if [ "${normauto}" != "chr1" ] && [ -z "${noCHR}" ]; then
	region1=$(grep $normauto $file.exons.bed | awk -F"\t" 'BEGIN{a=1000000000;b=0}{c=$1; if($2<a)a=$2; if($3>b)b=$3;} END{print c ":" a "-" b}')
fi
if [ "${normX}" != "chrX" ] && [ -z "${noCHR}" ]; then
	region2=$(grep $normX $file.exons.bed | awk -F"\t" 'BEGIN{a=1000000000;b=0}{c=$1; if($2<a)a=$2; if($3>b)b=$3;} END{print c ":" a "-" b}')
fi
if [ "${normauto}" != "chr1" ] && [ ! -z "${noCHR}" ]; then
	region1=$(grep $normauto $file.exons.bed | awk -F"\t" 'BEGIN{a=1000000000;b=0}{c=$1; if($2<a)a=$2; if($3>b)b=$3;} END{print c ":" a "-" b}' | sed 's/chr//g')
fi
if [ "${normX}" != "chrX" ] && [ ! -z "${noCHR}" ]; then
	region2=$(grep $normX $file.exons.bed | awk -F"\t" 'BEGIN{a=1000000000;b=0}{c=$1; if($2<a)a=$2; if($3>b)b=$3;} END{print c ":" a "-" b}' | sed 's/chr//g')
fi
rm -f $out/temp/region.grep.temp* $file.exons.*

# large regions for WGS
if [ ! -z "${WGS}" ]; then
	call=$script/data/hg19_regions-call.WGS.bed
fi
# small regions for WES
if [ -z "${WGS}" ]; then
	call=$script/data/hg19_regions-call.WES.bed
fi

# removing "chr" if --noCHR is used
if [ ! -z "${noCHR}" ]; then
    sed -i -e 's/chr//g' $out/00_raw-sequences/0.region.${normauto}.${normX}.bed
fi
if [ ! -z "${noCHR}" ] && [ ! -z "${hg38}" ]; then
    sed -i -e 's/chr//g' $out/00_raw-sequences/0.region.hg38.${normauto}.${normX}.bed
fi

if [ "$bam" == "No" ]
then
	echo "Step 2: Greping sequences in FASTQ files"

	IFS=$'\r\n' GLOBIGNORE='*' command eval 'patients=($(cut -f1 $input))'
	IFS=$'\r\n' GLOBIGNORE='*' command eval 'fastq1=($(cut -f2 $input))'
	IFS=$'\r\n' GLOBIGNORE='*' command eval 'fastq2=($(cut -f3 $input))'
	n=$((${#patients[@]}-1))

	for i in `seq 0 $n`
	do
		pat=${patients[$i]}

		echo "   Processing sample: $pat"
		# greping all sequences
		if [ ! -f $out/00_raw-sequences/$pat.sel.fastq.gz ]
		then
			# for paired-end
			if [ ! -z "${fastq2[$i]}" ]
			then
				# greping IDs above matching sequences
				zcat ${fastq1[$i]} | grep -B1 -Ff $out/00_raw-sequences/0.region.grep.${normauto}.${normX}.fa > $out/00_raw-sequences/$pat.read1.tsv &
				pid1=$!
				zcat ${fastq2[$i]} | grep -B1 -Ff $out/00_raw-sequences/0.region.grep.${normauto}.${normX}.fa > $out/00_raw-sequences/$pat.read2.tsv &
				pid2=$!
				wait $pid2
				wait $pid1

				cat $out/00_raw-sequences/$pat.read1.tsv $out/00_raw-sequences/$pat.read2.tsv | grep "@" | cut -f1 -d" " | sort | uniq > $out/00_raw-sequences/$pat.read.tsv

				zcat ${fastq1[$i]} | grep --no-group-separator -A3 -Ff $out/00_raw-sequences/$pat.read.tsv > $out/00_raw-sequences/$pat.1.fq &
				pid1=$!
				zcat ${fastq2[$i]} | grep --no-group-separator -A3 -Ff $out/00_raw-sequences/$pat.read.tsv > $out/00_raw-sequences/$pat.2.fq &
				pid2=$!
				wait $pid2
				wait $pid1

				$script/data/bbmap/bbmerge.sh in1=$out/00_raw-sequences/$pat.1.fq in2=$out/00_raw-sequences/$pat.2.fq out=$out/00_raw-sequences/$pat.merge.fq outu=$out/00_raw-sequences/$pat.nomerge.fq >/dev/null 2>&1
				cat $out/00_raw-sequences/$pat.merge.fq $out/00_raw-sequences/$pat.nomerge.fq > $out/00_raw-sequences/$pat.all.fq
				grep --no-group-separator -B1 -A2 -Ff $out/00_raw-sequences/0.region.grep.${normauto}.${normX}.fa $out/00_raw-sequences/$pat.all.fq > $out/00_raw-sequences/$pat.sel.fastq

			fi
			
			# for single-end
			if [ -z "${fastq2[$i]}" ]
			then
				zcat ${fastq1[$i]} | grep --no-group-separator -B1 -A2 -Ff $out/00_raw-sequences/0.region.grep.${normauto}.${normX}.fa > $out/00_raw-sequences/$pat.sel.fastq
			fi
			
			gzip $out/00_raw-sequences/$pat.sel.fastq
			rm -f $out/00_raw-sequences/${pat}.*.tsv $out/00_raw-sequences/${pat}.*.fq
		fi
	done
fi

if [ "$bam" == "Yes" ]
then
	echo "Step 2: Extracting reads from BAM files"

	IFS=$'\r\n' GLOBIGNORE='*' command eval 'patients=($(cut -f1 $input))'
	IFS=$'\r\n' GLOBIGNORE='*' command eval 'bams=($(cut -f2 $input))'
	n=$((${#patients[@]}-1))

	for i in `seq 0 $n`
	do
		pat=${patients[$i]}

		echo "   Processing sample: $pat"
		# greping all sequences
		if [ ! -f $out/00_raw-sequences/$pat.sel.fastq.gz ]
		then
			if [ ! -z "${hg38}" ]; then
				samtools view -h -M -L $out/00_raw-sequences/0.region.hg38.${normauto}.${normX}.bed ${bams[$i]} > $out/00_raw-sequences/${pat}.sam
			fi
			if [ -z "${hg38}" ]; then
				samtools view -h -M -L $out/00_raw-sequences/0.region.${normauto}.${normX}.bed ${bams[$i]} > $out/00_raw-sequences/${pat}.sam
			fi
			samtools view -bS $out/00_raw-sequences/${pat}.sam  > $out/00_raw-sequences/${pat}.bam
			samtools sort -o $out/00_raw-sequences/${pat}.sort.bam -n $out/00_raw-sequences/${pat}.bam 
			bedtools bamtofastq -i $out/00_raw-sequences/${pat}.sort.bam -fq $out/00_raw-sequences/${pat}.1.fq -fq2 $out/00_raw-sequences/${pat}.2.fq >/dev/null 2>&1

			bash $script/data/bbmap/bbmerge.sh in1=$out/00_raw-sequences/$pat.1.fq in2=$out/00_raw-sequences/$pat.2.fq out=$out/00_raw-sequences/$pat.merge.fq outu=$out/00_raw-sequences/$pat.nomerge.fq  >/dev/null 2>&1
			cat $out/00_raw-sequences/$pat.merge.fq $out/00_raw-sequences/$pat.nomerge.fq > $out/00_raw-sequences/$pat.all.fq
			grep --no-group-separator -B1 -A2 -Ff $out/00_raw-sequences/0.region.grep.${normauto}.${normX}.fa $out/00_raw-sequences/$pat.all.fq > $out/00_raw-sequences/$pat.sel.fastq
			gzip $out/00_raw-sequences/$pat.sel.fastq
			rm $out/00_raw-sequences/*am $out/00_raw-sequences/*fq
		fi
	done
fi


echo "Step 3: Determining coverage and haplotypes"
n=$((${#patients[@]}-1))
for i in `seq 0 $n`
do
	pat=${patients[$i]}
	echo "   Processing sample: $pat"

	file=$out/01_bam/$pat.bam
	
	if [ ! -f $file ]
	then
		# recreate BAM file: align selected reads and then convert SAM to BAM
		RG=$(echo "@RG\tID:$pat\tSM:$pat\tPL:Illumina")
		$bwa mem -M -C -R $RG $ref1 $out/00_raw-sequences/$pat.sel.fastq.gz > $out/temp/$pat.sam 2>> $out/08_logs/BWA-genome.txt
		$java -jar $picard SortSam -VALIDATION_STRINGENCY SILENT -SORT_ORDER coordinate -I $out/temp/$pat.sam -O $out/01_bam/$pat.bam -TMP_DIR $out/temp >> $out/08_logs/picard-genome.txt 2>&1
		$samtools index $out/01_bam/$pat.bam
		rm $out/temp/$pat.sam
	fi

	if [ ! -f $out/05_haplotypes/$pat.haplotypes.exon3.tsv ]
	then

		if [ ! -f $out/00_raw-sequences/$pat.seq.all-both.tsv ]
		then
			# extract sequences mapping to OPN1 region
			$samtools view $out/01_bam/$pat.bam "chrX:153390252-153510581" | cut -f10 > $out/temp/$pat.seq.all.tsv
			cat $out/temp/$pat.seq.all.tsv | tr ACGTacgt TGCAtgca | rev > $out/temp/$pat.seq.all-rev.tsv
			cat $out/temp/$pat.seq.all.tsv $out/temp/$pat.seq.all-rev.tsv > $out/temp/$pat.seq.all-both.tsv
			cp $out/temp/$pat.seq.all-both.tsv $out/00_raw-sequences
		fi

		# exon3 OPN1 sense based on two regions with changes in haplotypes
		grep -P "GTGAGATTTGATGCCAAGCTGGCCATC|TGCAAGCCCTTTGGCAA" $out/00_raw-sequences/$pat.seq.all-both.tsv > $out/temp/$pat.seq.hap-all.tsv
		# grep sequence that are long enough on both sides, first left and then right
		grep -P "GTGGCTGGTGGT|ATGGATGGTGGT|GTGGATGGTGGT|ATGGCTGGTGGT" $out/temp/$pat.seq.hap-all.tsv | grep -P "CTTCTCCTGGATCTGGT|CTTCTCCTGGGTCTGGT|CTTCTCCTGGGTCTGGG|CTTCTCCTGGATCTGGG" > $out/temp/$pat.sequences.tsv
		# R script for sequence alignement and haplotype calling
		if [ -f $out/temp/$pat.sequences.tsv ] && [ -s $out/temp/$pat.sequences.tsv ]
		then
			$Rscript $script/sequence-analysis_exon3.R $out/temp/$pat.sequences.tsv $out/05_haplotypes/$pat.haplotypes.exon3.tsv $out/05_haplotypes/$pat.haplotypesDNA.exon3.tsv
		fi
		# remove temporary files
		rm $out/temp/$pat.seq.hap*

		# exon4 haplotypes
		grep -P "CATGGTCACCTGCTGCATCA|TGCTCTGCTACCTCCAAGTG" $out/00_raw-sequences/$pat.seq.all-both.tsv > $out/temp/$pat.seq.hap-all.tsv
		grep -P "CTGCTGCATCA" $out/temp/$pat.seq.hap-all.tsv | grep -P "TGCTCTGCTAC" > $out/temp/$pat.sequences.tsv
		if [ -f $out/temp/$pat.sequences.tsv ] && [ -s $out/temp/$pat.sequences.tsv ]
		then
			$Rscript $script/sequence-analysis_exon4.R $out/temp/$pat.sequences.tsv $out/05_haplotypes/$pat.haplotypes.exon4.tsv
		fi
		rm $out/temp/$pat.seq.hap*

		# exon2 haplotypes
		grep -P "TGCATCCGTCTTCACAAATGGGC" $out/00_raw-sequences/$pat.seq.all-both.tsv > $out/temp/$pat.seq.hap-all.tsv
		grep -P "TGCATCCGTCTTCACAAATGGGCT" $out/temp/$pat.seq.hap-all.tsv | grep -P "TTGTGAACCAGGTCT" > $out/temp/$pat.sequences.tsv
		if [ -f $out/temp/$pat.sequences.tsv ] && [ -s $out/temp/$pat.sequences.tsv ]
		then
			$Rscript $script/sequence-analysis_exon2.R $out/temp/$pat.sequences.tsv $out/05_haplotypes/$pat.haplotypes.exon2.tsv
		fi
		rm $out/temp/$pat.seq.hap*

		# exon5 haplotypes
		grep -P "TCGCATGCTTTGCTGCTGCCAAC" $out/00_raw-sequences/$pat.seq.all-both.tsv > $out/temp/$pat.seq.hap-all.tsv
		grep -P "TCTGCTGGGGACC" $out/temp/$pat.seq.hap-all.tsv | grep -P "CTGCCCTGCCGGCCT" > $out/temp/$pat.sequences.tsv
		if [ -f $out/temp/$pat.sequences.tsv ] && [ -s $out/temp/$pat.sequences.tsv ]
		then
			$Rscript $script/sequence-analysis_exon5.R $out/temp/$pat.sequences.tsv $out/05_haplotypes/$pat.haplotypes.exon5.tsv
		fi
		rm $out/temp/$pat.seq.hap*

	fi

	if [ ! -f $out/04_coverage/$pat.coverage.tsv ]
	then

		if [ ! -f $out/00_raw-sequences/$pat.seq.all-both.tsv ]
		then
			# extract sequences mapping to OPN1 region
			$samtools view $out/01_bam/$pat.bam "chrX:153390252-153510581" | cut -f10 > $out/temp/$pat.seq.all.tsv
			cat $out/temp/$pat.seq.all.tsv | tr ACGTacgt TGCAtgca | rev > $out/temp/$pat.seq.all-rev.tsv
			cat $out/temp/$pat.seq.all.tsv $out/temp/$pat.seq.all-rev.tsv > $out/temp/$pat.seq.all-both.tsv
			cp $out/temp/$pat.seq.all-both.tsv $out/00_raw-sequences
		fi

		# look at coverage
		grep CACAGGCCAGTATAAAGCGCCGTGAC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon1-part1" "\t" $1}' > $out/04_coverage/$pat.coverage.tsv
		grep CACTGGCCGGTATAAAGCACCGTGAC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon1-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CTCAGGTGATGCGCCAGGGCCGGCT $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon1-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CTCAGGTGACGCACCAGGGCCGGCT $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon1-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep GATGATCTTTGTGGTCACTGCATCCGTC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon2-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep GATGATCTTTGTGGTCATTGCATCCGTC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon2-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep GTGAACCAGGTCTCTGGCTACTTCGTGCTGGG $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon2-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep GTGAACCAGGTCTATGGCTACTTCGTGCTGGG $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon2-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep GTGAGATTTGATGCCAAGCTGGCCATC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1-exon3-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep TGTGTGGACAGCCCCGCCCATCTTTGGTTGGAGCAG $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1-exon3-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep CACCTGCTGCATCATCCCACTCG $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon4-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CACCTGCTGCATCACCCCACTCA $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon4-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CTATCATCATGCTCTGCTACCTC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon4-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep GCATCATCGTGCTCTGCTACCTC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MWMW2-exon4-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep TACTGCGTCTGCTGGGGACCC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon5-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep TTCTGCTTCTGCTGGGGACCA $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MW-exon5-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep TTCTGCTTCTGCTGGGGACCC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MW2-exon5-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep CTACACCTTCTTCGCATGCTT $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1LW-exon5-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep TACGCCTTCTTCGCATGCTTT $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MW-exon5-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep TACGCCTTCTTCGCATGCTTT $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1MW2-exon5-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		grep GAAACTGCATCTTGCAGCTTT $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1-exon6-part1" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		grep AGGTCTCATCTGTGTCCTCGG $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" '{print "OPN1-exon6-part2" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

		if [ "${normX}" == "chrX" ] && [ -z "${noCHR}" ]; then
			$samtools idxstats ${bams[$i]} | cut -f1,3 | grep -P "chrX\t" | cut -f2 |awk -F"\t" '{print "Control-chrX" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		fi
		if [ "${normauto}" == "chr1" ] && [ -z "${noCHR}" ]; then
			$samtools idxstats ${bams[$i]} | cut -f1,3 | grep -P "chr1\t" | cut -f2 | awk -F"\t" '{print "Control-autosome" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		fi
		if [ "${normX}" == "chrX" ] && [ ! -z "${noCHR}" ]; then
			$samtools idxstats ${bams[$i]} | cut -f1,3 | grep -P "X\t" | cut -f2 |awk -F"\t" '{print "Control-chrX" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		fi
		if [ "${normauto}" == "chr1" ] && [ ! -z "${noCHR}" ]; then
			$samtools idxstats ${bams[$i]} | cut -f1,3 | grep -P "^1\t" | cut -f2 | awk -F"\t" '{print "Control-autosome" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		fi

		if [ "${normX}" != "chrX" ]; then
			$samtools view $out/01_bam/$pat.bam $region2 | wc -l | awk -F"\t" '{print "Control-chrX" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		fi
		if [ "${normauto}" != "chr1" ]; then
			$samtools view $out/01_bam/$pat.bam $region1 | wc -l | awk -F"\t" '{print "Control-autosome" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		fi

		$samtools view $out/01_bam/$pat.bam $region3 | wc -l | awk -F"\t" '{print "LCR" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv
		$samtools view $out/01_bam/$pat.bam "chrX:153405983-153406393" | wc -l | awk -F"\t" '{print "LCR-small" "\t" $1}' >> $out/04_coverage/$pat.coverage.tsv

	fi
	if [ ! -f $out/04_coverage/$pat.coverage.mut.tsv ]
	then

		if [ ! -f $out/00_raw-sequences/$pat.seq.all-both.tsv ]
		then
			# extract sequences mapping to OPN1 region
			$samtools view $out/01_bam/$pat.bam "chrX:153390252-153510581" | cut -f10 > $out/temp/$pat.seq.all.tsv
			cat $out/temp/$pat.seq.all.tsv | tr ACGTacgt TGCAtgca | rev > $out/temp/$pat.seq.all-rev.tsv
			cat $out/temp/$pat.seq.all.tsv $out/temp/$pat.seq.all-rev.tsv > $out/temp/$pat.seq.all-both.tsv
			cp $out/temp/$pat.seq.all-both.tsv $out/00_raw-sequences
		fi
		
		# mutations in exon 4
		grep CGCGGCCCAGACGTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAT $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "C203R-LW" "\t" $1}' > $out/04_coverage/$pat.coverage.mut.tsv
		grep CGCGGCCCAGACGTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "C203R-MW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		grep TGCGGCCCAGACGTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAT $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "C203-LW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
		grep TGCGGCCCAGACGTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCAC $out/00_raw-sequences/$pat.seq.all-both.tsv | wc -l | awk -F"\t" -v pat="$pat" '{print pat "\t" "C203-MW" "\t" $1}' >> $out/04_coverage/$pat.coverage.mut.tsv
	fi

	rm -f $out/temp/$pat.* $out/00_raw-sequences/$pat.seq.all-both.tsv
	
done

echo "Step 4: Analysis of rare variants"
for pat in "${patients[@]}"
do
	echo "   Processing sample: $pat"
	if [ ! -f $out/03_variants/$pat.variants.tsv ] || [ ! -f $out/01_bam/$pat.OPN1LWMW.bam ] || [ ! -f $out/02_vcf/$pat.norm.vcf ]
	then
		# take reads from BAM (OPN1 region) and create FASTQs
		$samtools view -b -h $out/01_bam/$pat.bam "chrX:153390252-153510581" > $out/temp/$pat.region.bam
		$samtools index $out/temp/$pat.region.bam
		$samtools sort $out/temp/$pat.region.bam > $out/temp/$pat.region.sorted.bam
		$samtools index $out/temp/$pat.region.sorted.bam
		$bedtools bamtofastq -i $out/temp/$pat.region.sorted.bam -fq $out/temp/$pat.region.fastq

		# alignement, sam to bam and indexing
		RG=$(echo "@RG\tID:$pat\tSM:$pat\tPL:Illumina")
		$bwa mem -M -C -R $RG $ref2 $out/temp/$pat.region.fastq > $out/temp/$pat.region.sam 2>> $out/08_logs/BWA-region.txt
		$java -jar $picard SortSam -VALIDATION_STRINGENCY SILENT -SORT_ORDER coordinate -I $out/temp/$pat.region.sam -O $out/temp/$pat.OPN1LWMW.bam -TMP_DIR $out/temp 2>> $out/08_logs/picard-region.txt
		$samtools index $out/temp/$pat.OPN1LWMW.bam
		cp $out/temp/$pat.OPN1LWMW.bam* $out/01_bam/

		# variant calling
		$java -jar $gatk HaplotypeCaller \
		-I $out/temp/$pat.OPN1LWMW.bam \
		-R $ref2 \
		-O $out/temp/$pat.vcf \
		-G StandardAnnotation \
		-G StandardHCAnnotation \
		-L $call 2>> $out/08_logs/HaplotypeCaller.txt

		# VCF processing
		$bcftools norm -m -both --fasta-ref $ref2 -o $out/02_vcf/$pat.norm.vcf $out/temp/$pat.vcf 2>> $out/08_logs/bcftools.txt
		grep -P "chrX\t" $out/02_vcf/$pat.norm.vcf | awk -F"\t" -v pat="$pat" '{split($10,a,":"); split(a[2],b,","); print pat "\t" $1 "\t" $2 "\t" $4 "\t" $5 "\t" b[1] "\t" b[2] "\t" $1 "-" $2 "-" $4 "-" $5}' > $out/03_variants/$pat.variants.tsv
	fi
done

echo "Step 5: Batch processing"

# merge variants in one file
echo -e "ID\tchr\tpos\tref\talt\treads-ref\treads-alt\tchange" > $out/03_variants/0.variants-all.tsv
cat $out/03_variants/*.variants.tsv >> $out/03_variants/0.variants-all.tsv

echo -e "ID\tchr\tpos\tref\talt\treads-ref\treads-alt\tchange" > $out/03_variants/0.variants-rare.tsv
cat $out/03_variants/*.variants.tsv | grep -v -Ff $script/data/gnomAD-frequent.tsv >> $out/03_variants/0.variants-rare.tsv

echo -e "ID\tType\tReads" > $out/04_coverage/0.coverage.mut.ALL.tsv
cat $out/04_coverage/*.coverage.mut.tsv >> $out/04_coverage/0.coverage.mut.ALL.tsv

# extract pathogenic and rare haplotypes for exon 3
echo -e "ID\tHaplotype\tReads\tTotal\tPerc" > $out/05_haplotypes/0.exon3.PROT.all.tsv
grep -P "A|I|V|S" $out/05_haplotypes/*.haplotypes.exon3.tsv | awk -F"\t" '{split($1,a,":"); n=split(a[1],b,"/"); split(b[n],c,"."); print c[1] "\t" a[2] "\t" $2 "\t" $3 "\t" $4}' >> $out/05_haplotypes/0.exon3.PROT.all.tsv

# extract rare DNA haplotypes for exon 3
echo -e "ID\tHaplotype\tReads\tTotal\tPerc" > $out/05_haplotypes/0.exon3.DNA.all.tsv
grep -v "Haplotype" $out/05_haplotypes/*.haplotypesDNA.exon3.tsv | awk -F"\t" '{split($1,a,":"); n=split(a[1],b,"/"); split(b[n],c,"."); print c[1] "\t" a[2] "\t" $2 "\t" $3 "\t" $4}' >> $out/05_haplotypes/0.exon3.DNA.all.tsv

# extract pathogenic haplotypes based on splicing assay
awk -F'\t' -v OFS='\t' 'NR==FNR {map[$2] = $1; next }{ print $0, ( $2 in map ? map[$2] : "Haplotype" )}' $script/haplo_exon3_splicing.tsv $out/05_haplotypes/0.exon3.DNA.all.tsv > $out/05_haplotypes/temp.exon3.tsv
awk -F"\t" '{if($3>=0.70 && $3<2) print $2}' $script/haplo_exon3_splicing.tsv > $out/05_haplotypes/temp.spli.tsv
echo -e "Haplotype" >> $out/05_haplotypes/temp.spli.tsv
echo -e "ID\tHaplotypeDNA\tHaplotype\tReads\tTotal\tPerc" > $out/05_haplotypes/0.exon3.pathogenic.tsv
grep -Ff $out/05_haplotypes/temp.spli.tsv $out/05_haplotypes/temp.exon3.tsv | awk -F"\t" '{print $1 "\t" $2 "\t" $6 "\t" $3 "\t" $4 "\t" $5}' | grep -v "Haplotype" >> $out/05_haplotypes/0.exon3.pathogenic.tsv
cut -f1,3,4,5,6 $out/05_haplotypes/0.exon3.pathogenic.tsv > $out/05_haplotypes/0.exon3.PROT.pathogenic.tsv
cut -f1,2,4,5,6 $out/05_haplotypes/0.exon3.pathogenic.tsv > $out/05_haplotypes/0.exon3.DNA.pathogenic.tsv
file="$out/05_haplotypes/0.exon3.PROT.pathogenic.tsv"
if [ "$(wc -l < "$file")" -lt 2 ]; then
    echo -e "NA\tNA\tNA\tNA\tNA" >> "$file"
fi
file="$out/05_haplotypes/0.exon3.DNA.pathogenic.tsv"
if [ "$(wc -l < "$file")" -lt 2 ]; then
    echo -e "NA\tNA\tNA\tNA\tNA" >> "$file"
fi
file="$out/05_haplotypes/0.exon3.pathogenic.tsv"
if [ "$(wc -l < "$file")" -lt 2 ]; then
    echo -e "NA\tNA\tNA\tNA\tNA\tNA" >> "$file"
fi
rm $out/05_haplotypes/temp.exon3.tsv $out/05_haplotypes/temp.spli.tsv

# extract other exon haplotypes
echo -e "ID\tHaplotype\tExon\tReads\tTotal\tPerc" > $out/05_haplotypes/0.exon2.PROT.all.tsv
echo -e "ID\tHaplotype\tExon\tReads\tTotal\tPerc" > $out/05_haplotypes/0.exon4.PROT.all.tsv
echo -e "ID\tHaplotype\tExon\tReads\tTotal\tPerc" > $out/05_haplotypes/0.exon5.PROT.all.tsv
grep -P -v "Haplotype" $out/05_haplotypes/*.haplotypes.exon2.tsv | awk -F"\t" '{split($1,a,":"); n=split(a[1],b,"/"); split(b[n],c,"."); print c[1] "\texon2\t" a[2] "\t" $2 "\t" $3 "\t" $4}' >> $out/05_haplotypes/0.exon2.PROT.all.tsv
grep -P -v "Haplotype" $out/05_haplotypes/*.haplotypes.exon4.tsv | awk -F"\t" '{split($1,a,":"); n=split(a[1],b,"/"); split(b[n],c,"."); print c[1] "\texon4\t" a[2] "\t" $2 "\t" $3 "\t" $4}' >> $out/05_haplotypes/0.exon4.PROT.all.tsv
grep -P -v "Haplotype" $out/05_haplotypes/*.haplotypes.exon5.tsv | awk -F"\t" '{split($1,a,":"); n=split(a[1],b,"/"); split(b[n],c,"."); print c[1] "\texon5\t" a[2] "\t" $2 "\t" $3 "\t" $4}' >> $out/05_haplotypes/0.exon5.PROT.all.tsv
echo -e "ID\tExon\tHaplotype\tReads\tTotal\tPerc" > $out/05_haplotypes/0.exon2-4-5.PROT.rare.tsv
cat $out/05_haplotypes/0.exon2.PROT.all.tsv $out/05_haplotypes/0.exon4.PROT.all.tsv $out/05_haplotypes/0.exon5.PROT.all.tsv | grep -v -P "Haplotype|\tIVY\t|\tTIS\t|\tTSV\t|\tIAM\t|\tIFYVTAY\t|\tVLFFAPF\t" | awk -F"\t" '{if($4>4) print $0}' >> $out/05_haplotypes/0.exon2-4-5.PROT.rare.tsv

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

echo "Step 6: Creating plots and events files"
$Rscript $script/analysis_plots.R $out/04_coverage/0.coverage-ALL.tsv $out/05_haplotypes/0.exon3.PROT.all.tsv $out/05_haplotypes/0.exon3.pathogenic.tsv $out/05_haplotypes/0.exon3.DNA.all.tsv $out/03_variants/0.variants-rare.tsv $script/data/HGVS-gnomAD.RData $out $out/04_coverage/0.coverage.mut.ALL.tsv $out/03_variants/0.variants-all.tsv $malesonly $script/data/Opsin-logo.png $sexfile

rm -rf $out/temp $out/03_variants/0.variants-rare.tsv $out/03_variants/0.variants-all.tsv

