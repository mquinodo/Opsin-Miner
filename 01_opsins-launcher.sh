
### Made by Mathieu Quinodoz
### October 2022

here=/home/mquinodo
script=$here/SYNO/scripts_NGS_analysis/Scripts_WES_analysis/opsins/v5
NAS=$here/SYNO/WES/EXOMES

batch=Carmen_opsin

out=$NAS/$batch/11_opsins-v3
mkdir -p $out

input=$out/$batch.input.tsv
a=$(ls $NAS/$batch/*/original_data/*_1.fastq.gz | cut -d"/" -f8 | tr '\n' ',')
names=${a::-1}
IFS=',' read -r -a patients <<< "$names"
rm -f $input
touch $input
for pat in "${patients[@]}"
do
	printf "$pat\t" >> $input
	printf "$NAS/$batch/$pat/original_data/${pat}_1.fastq.gz\t" >> $input
	printf "$NAS/$batch/$pat/original_data/${pat}_2.fastq.gz\n" >> $input
done

picard=/usr/local/bin/picard.jar
gatk=/usr/local/bin/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar

bash $script/02_opsins-main.sh --script $script --NAS $NAS --out $out --input $input --picard $picard --gatk $gatk #>> $out/log.txt 2>&1

######

here=/home/mquinodo/SYNO/WES/EXOMES

ls -lrt $here/*/11_opsins-v3/0.e*

cat $here/*/11_opsins-v3/0.e* | grep -v -P "T1527A|MAPK1151|CRMM1|CRMM2|CM2005|\-B|\-T|PK|KN-" | sort | uniq > $here/all.tsv
grep -P "#|High" $here/all.tsv > $here/all-HQ.tsv
cat $here/*/11_opsins-v3/0.s* | grep -v -P "T1527A|MAPK1151|CRMM1|CRMM2|CM2005|\-B|\-T|PK|KN-" > $here/all-sex.tsv

cat $here/all-sex.tsv | grep Male | wc -l
cat $here/all-HQ.tsv | grep Male | grep deletion | grep -v Event | cut -f1 | sort | uniq | wc -l

cat $here/all-sex.tsv | grep Female | wc -l
cat $here/all-HQ.tsv | grep Female | grep deletion | grep -v Event | cut -f1 | sort | uniq | wc -l


#######


### Made by Mathieu Quinodoz
### October 2022

here=/home/mquinodo
script=$here/SYNO/scripts_NGS_analysis/Scripts_WES_analysis/opsins/v5
NAS=$here/SYNO/WES/EXOMES

for batch in Carmen_Igor Carmen_opsin CeGat_2020-03 CeGat_2020-04 CeGat_2020-07 CeGat_2020-08 CeGat_2020-10 CeGat_2020-11 CeGat_2020-12 CeGat_2021-01 CeGat_2021-02 CeGat_2021-03 CeGat_2021-05 CeGat_2021-06 CeGat_2021-07 CeGat_2021-10 CeGat_2021-11 CeGat_2021-12 CeGat_2022-02 CeGat_2022-03 CeGat_2022-04 CeGat_2022-07 CeGat_2022-09 ESTONIA ESTONIA_2 ESTONIA_3 IRD_GREEK KOSTA_CRD LISBON NEMIPS NOVOGENE_180913 NOVOGENE_181113 Novogene_Jan2019 Novogene_Jan2020 Novogene_May2019 Novogene_Oct2019 Novogene_Sept2019
do

	out=$NAS/$batch/11_opsins-v3
	mkdir -p $out

	input=$out/$batch.input.tsv

	picard=/usr/local/bin/picard.jar
	gatk=/usr/local/bin/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar

	bash $script/02_opsins-main.sh --script $script --NAS $NAS --out $out --input $input --picard $picard --gatk $gatk #>> $out/log.txt 2>&1

done



