#!/bin/bash

# $1 - sample binning dir, $2 - sample taira file

cd $1

mkdir -p bins
mkdir -p ebins
mkdir -p stats

for f in $(ls metabat); do
    cp metabat/$f bins/${f%.fa}.fna
done

for f in $(ls vamb/bins); do
    cp vamb/bins/$f bins
done

for f in $(ls semibin/output_bins); do
    cp semibin/output_bins/$f bins/${f%.fa.gz}.fna.gz
done

gzip -d bins/*.gz

seqkit stats bins/*.fna -o $1/stats/seqkit.stats.tsv

for f in $(awk '{print $1}' stats/seqkit.stats.tsv | tail -n +2);do 
    grep -wFf <(grep ">" $f | tr -d '>') $2 | \
        cut -f 2 | sort | uniq -c | sort -nr -k 1 | head -n 1
done | paste <(tail -n +2 stats/seqkit.stats.tsv) - > stats/seqkit.stats2.tsv

for f in $(grep "eukarya" stats/seqkit.stats2.tsv | tr -d ',' | awk '{if ($5 > 5000000) print $1}'); do
    mv $f ebins
done