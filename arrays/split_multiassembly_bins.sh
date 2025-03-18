#!/bin/bash

FILE=$1
BASE=${FILE%.fna}

# split multiassembly into separate sample files
for s in `grep -oE ">[0-9]+" $1 | tr -d '>' | sort | uniq`; do
    seqkit grep -rp "$s*" $FILE > ${BASE}_$s.fna
done

# remove samples with less than 5M nucleotides
for f in `wc -m ${BASE}*.fna | awk '$1 < 5000000' | awk '{ print $2 }'`; do
    rm -f $f
done


