#!/bin/bash

PICARD_PATH=/mnt/lustre/home/shyamg/tools/Picard/picard-tools-1.54
JAVA=/usr/bin/java
SAMTOOLS=/data/tools/samtools_new/samtools/samtools

# read raw bam file from command line
BAMFILE=$1
LOGPATH=`dirname $BAMFILE`
TAG=`basename $BAMFILE`
tmpdir=${LOGPATH}/tmp
LOGPATH=${LOGPATH}/${TAG}

# convert to bam, sort and index
echo "Sorting aligned reads..."
SORTED_BAMFILE=`echo $BAMFILE | sed 's/.bam/.sort.bam/'`
$JAVA -jar -Djava.io.tmpdir=$tmpdir ${PICARD_PATH}/SortSam.jar SO=coordinate I=${BAMFILE} O=${SORTED_BAMFILE} 
echo "Indexing aligned reads..."
$SAMTOOLS index ${SORTED_BAMFILE}

remove mitochondria
BAMFILE=`echo $SORTED_BAMFILE | sed 's/.bam/.nuc.bam/'`
MITOFILE=`echo $SORTED_BAMFILE | sed 's/.bam/.chrM.bam/'`
chrs=`samtools view -H ${SORTED_BAMFILE} | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM`
echo "Removing mitochondial reads ..."
samtools view -b $SORTED_BAMFILE `echo $chrs` > $BAMFILE
samtools view -b $SORTED_BAMFILE chrM > $MITOFILE

# samtools rmdup
echo "Removing duplicates..."
RMDUP_BAMFILE=`echo $BAMFILE | sed 's/.bam/.rmdup.bam/'`
$JAVA -jar -Djava.io.tmpdir=$tmpdir ${PICARD_PATH}/MarkDuplicates.jar INPUT=$BAMFILE OUTPUT=$RMDUP_BAMFILE METRICS_FILE=${LOGPATH}.dups.log REMOVE_DUPLICATES=true

# index
echo "Indexing reads after removing duplicates ..."
$SAMTOOLS index $RMDUP_BAMFILE

# histogram file
echo "Generating fragment size distribution for QC ..."
$JAVA -jar -Djava.io.tmpdir=$tmpdir ${PICARD_PATH}/CollectInsertSizeMetrics.jar I=$RMDUP_BAMFILE O=${LOGPATH}.hist_data.log H=${LOGPATH}.hist_graph.pdf W=1000 STOP_AFTER=50000000
#rm -rf $tmpdir

$JAVA -jar -Djava.io.tmpdir=$tmpdir ${PICARD_PATH}/CollectInsertSizeMetrics.jar I=$MITOFILE O=${LOGPATH}.${MITOFILE}.hist_data.log H=${LOGPATH}.${MITOFILE}.hist_graph.pdf W=1000 STOP_AFTER=50000000
rm -rf $tmpdir
