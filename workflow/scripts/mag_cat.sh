#!/bin/bash

SAMPLE=$(basename $1)

for f in `ls $1/bins_das/das_DASTool_bins`; do
    sed -i "s|>|>${SAMPLE}_|g" $f
done

for f in `awk '$10 > 49 && $11 < 10' $1/bins_das/das_DASTool_summary.tsv | cut -f 1`; do
    cp $1/bins_das/das_DASTool_bins/$f /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/prok/${SAMPLE}_$f
done


# This will strip fasta headers leaving just the ID and remove any comments. 
# Be careful because you must have uppercase nucleotide seqs and it could clobber other data!
# sed -i -e 's|\s[a-z0-9.=]*||g' $FILE