#!/usr/bin/env bash

# Run mmseqs2 using 

# Parameters:

# the sample id - 33000xxxxx
SAMPLE_ID=$1

# output directory (full path)
OUT=$2

# a temporary directory (full path)
TMP=$3

# mmseqs target DB 
TARGET_DB=$4

# n. threads
CPUS_PER_TASK=$5

# folder containing all the zipped sample annotations 
# e.g. /path/to/33000xxxxx.tar.gz
IMG_FOLDER=$6


if [ ! -d $TMP ]; then
    mkdir $TMP
fi

if [ ! -d $OUT ]; then
    mkdir $OUT
fi

zcat ${IMG_FOLDER}/${SAMPLE_ID}.tar.gz | tar xfO - ${SAMPLE_ID}.a.faa > $TMP/${SAMPLE_ID}.a.faa

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

rm -rf $TMP
rm -f $OUT/${SAMPLE_ID}.faa.mms2.[0-9]*
