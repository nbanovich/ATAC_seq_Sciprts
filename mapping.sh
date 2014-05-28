#!/bin/bash

fastq_R1=$1
genome_index=$2
fname=$(echo $fastq_R1 | sed 's/\_R1_001.fastq.gz//g')
fname_trim=$(echo $fname | sed 's/_L001//g' | sed 's/_L002//g' |sed 's/\//\t/g' | cut -f3 )
fname_name=$(echo $fname | sed 's/\//\t/g' | cut -f3 )
mkdir ./${fname_trim}
fastq_R2="${fname}_R2_001.fastq.gz" 
sam_name="./${fname_trim}/${fname_name}.sam"
bam_name="./${fname_trim}/${fname_name}.bam"
#sort_name="./${fname}/${fname_trim}.sorted"
/mnt/lustre/home/siddisis/bin/bowtie2 -X2000 -x $genome_index -1 $fastq_R1 -2 $fastq_R2 -S $sam_name 
samtools view -bS $sam_name > $bam_name
#samtools sort $bam_name $sort_name
exit
