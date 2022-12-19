
### Made by Mathieu Quinodoz
### October 2022

here=/home/mquinodo
script=$here/SYNO/scripts_NGS_analysis/Scripts_WES_analysis/opsins/Opsin-Miner-main
NAS=$here/SYNO/WES/EXOMES

batch=Carmen_opsin

out=$NAS/$batch/11_opsins-v5
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

bash $script/02_opsins-main.sh --script $script --out $out --input $input --picard $picard --gatk $gatk >> $out/log.txt 2>&1


