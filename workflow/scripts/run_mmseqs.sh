#!/usr/bin/env bash

# Run mmseqs2 using 

SAMPLE_ID=$1
OUT=$2
TMP=$3
TARGET_DB=$4
CPUS_PER_TASK=$5
IMG_FOLDER=$6


if [ ! -d $TMP ]; then
    mkdir $TMP
fi

zcat  ${IMG_FOLDER}/${SAMPLE_ID}.tar.gz | tar xfO - ${SAMPLE_ID}.a.faa > $TMP/${SAMPLE_ID}.a.faa

if [ ! -d $OUT ]; then
    mkdir $OUT
fi

QUERY=$TMP/${SAMPLE_ID}_DB

# Index query sequence file and create DB
mmseqs createdb $TMP/${SAMPLE_ID}.a.faa $QUERY


### Run mmseqs

# Run lowest common ancestor assignment
mmseqs taxonomy $QUERY ${TARGET_DB} $OUT/${SAMPLE_ID}.faa.mms2 $TMP \
  --lca-mode 3 \
  -s 4.0 \
  --threads ${CPUS_PER_TASK} \
  --blacklist 12908,28384 \
  --tax-lineage 1 \
  --compressed 1 \
  --remove-tmp-files 1

# Convert results to tsv
mmseqs createtsv \
  $QUERY ${OUT}/${SAMPLE_ID}.faa.mms2 \
  ${OUT}/${SAMPLE_ID}.faa.mms2.tsv


### Delete query DB and intermediate files

rm -f $TMP/${SAMPLE_ID}.a.faa
rm -f $QUERY*
rm -f $OUT/${SAMPLE_ID}.faa.mms2.[0-9]*
