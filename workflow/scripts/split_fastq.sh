#!/usr/bin/env bash

BAM=$1
EUKS=$2
PAIRED=$3
SINGLETONS=$4
THREADS=$5


#samtools index $BAM
LEN=$(samtools view -h $BAM | head -n 2 | tail -n 1 | cut -f 3 -d ":")

samtools view -h -F 12 -@ $THREADS \
    --regions-file <(sed "s|$|\t0\t$LEN|g" $EUKS) $BAM | \
samtools fastq -@ $THREADS -1 $PAIRED -2 $PAIRED -s $SINGLETONS
#rm -f $BAM.bai