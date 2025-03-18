#!/bin/bash

SAMPLE=$(basename $1)

#for f in `ls $1/bins_das/das_DASTool_bins`; do
#    sed -i "s|>|>${SAMPLE}_|g" $f
#done

for f in `awk '$10 > 49 && $11 < 10' $1/bins_das/das_DASTool_summary.tsv | cut -f 1 | grep "33000"`; do
    sed "s|>|>${SAMPLE}_|g" $1/bins_das/das_DASTool_bins/${f}.fa > /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/prok/${f}.fna
done

for f in `awk '$10 > 49 && $11 < 10' $1/bins_das/das_DASTool_summary.tsv | cut -f 1 | grep -v "33000"`; do
    sed "s|>|>${SAMPLE}_|g" $1/bins_das/das_DASTool_bins/${f}.fa > /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/prok/${SAMPLE}_${f}.fna
done

