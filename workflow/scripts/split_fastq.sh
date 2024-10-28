#!/usr/bin/env bash

BAM=$1
EUKS=$2
MAP=$3
PAIRED=$4
SINGLETONS=$5
THREADS=$6


samtools index $BAM
LEN=$(samtools view -h $BAM | head -n 2 | tail -n 1 | cut -f 3 -d ":")

samtools view -h -F 12 -@ $THREADS \
    --regions-file <(grep -wFf $EUKS $MAP | cut -f 2 | sed "s|$|\t0\t$LEN|g") $BAM | \
    samtools fastq -@ $THREADS -1 $PAIRED -2 $PAIRED -s $SINGLETONS
rm -f $BAM.bai