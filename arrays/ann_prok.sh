#!/bin/bash

source ~/start_conda

conda activate bio


cd /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue

for i in `seq 1 9987`; do 
    echo "$i/9987"
    f=$(ls prok/ | head -n $i | tail -n 1)
    grep -f \
        <(cat prok/$f | grep ">" | sed 's|_|\t|1' | cut -f 2) \
        /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/${f%%_*}/${f%%_*}.a.pfam.txt > \
        annotations/prok/pfam/${f%.fna}.pfam.txt

done

for i in `seq 1 9987`; do 
    echo "$i/9987"
    f=$(ls prok/ | head -n $i | tail -n 1)
    grep -f \
        <(cat prok/$f | grep ">" | sed 's|_|\t|1' | cut -f 2) \
        /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/${f%%_*}/${f%%_*}.a.gff > \
        annotations/prok/gff/${f%.fna}.gff
done

for i in `seq 1 9987`; do 
    echo "$i/9987"
    f=$(ls prok/ | head -n $i | tail -n 1)
    grep -f \
        <(cat prok/$f | grep ">" | sed 's|_|\t|1' | cut -f 2) \
        /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/${f%%_*}/${f%%_*}.a.gff > \
        annotations/prok/gff/${f%.fna}.gff
done

ls prok | \
    parallel -j 16 "f={};seqkit grep -rf <(grep '>' prok/{} | sed 's|_|\t|1' | cut -f 2 | sed 's|$|_*|g') /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/${f%%_*}/${f%%_*}.a.faa > annotations/prok/faa/\$f{%.*}.faa"


cut -f 1 prok_gtdb_tk_out/gtdbtk.bac120.summary.tsv | tail -n +2 | \
    parallel -j 16 "barrnap --kingdom bac --quiet prok/{}.fna | tail -n +2 >> annotations/prok/gff/{}.gff"

cut -f 1 prok_gtdb_tk_out/gtdbtk.arc83.summary.tsv | tail -n +2 | \
    parallel -j 16 "barrnap --kingdom arc --quiet prok/{}.fna | tail -n +2 >> annotations/prok/gff/{}.gff"

ls annotations/vir/gff | head -n 2 | parallel -j 16 "f={};for s in \$(sed 's|_|\t|1' annotations/vir/gff/{} | cut -f 1 | sort | uniq); do seqkit grep -f <(cut -f 9 annotations/vir/gff/{} | grep -oE 'ID=Ga[0-9_]*' | sed 's|ID=||g') /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/\$s/\$s.a.faa; done >> annotations/vir/faa/\${f%.*}.faa"

