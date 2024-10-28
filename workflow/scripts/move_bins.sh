#!/bin/bash

# $1 - sample binning dir, $2 - sample taira file

cd $1
SEQID=$(basename $1)
mkdir -p bins
mkdir -p ebins
mkdir -p bins_checkm
mkdir -p ebins_eukcc

gzip -d semibin/output_bins/*.gz

for f in $(ls metabat); do
    mv metabat/$f bins/${f%.fa}.fna
done

for f in $(ls vamb/bins); do
    mv vamb/bins/$f bins
done

for f in $(ls semibin/output_bins); do
    mv semibin/output_bins/$f bins/${f%.fa}.fna
done

seqkit stats bins/*.fna -o $1/bins_checkm/seqkit.stats.tsv

for f in $(awk '{print $1}' bins_checkm/seqkit.stats.tsv);do 
    grep -wFf <(grep ">" $f | tr -d '>') $2 | \
        cut -f 2 | sort | uniq -c | sort -nr -k 1 | head -n 1
done | paste <(tail -n +2 bins_checkm/seqkit.stats.tsv) - > bins_checkm/seqkit.stats2.tsv

for f in $(grep "eukarya" bins_checkm/seqkit.stats2.tsv | tr -d ',' | awk '{if ($5 > 5000000) print $1}'); do
    mv $f ebins
done