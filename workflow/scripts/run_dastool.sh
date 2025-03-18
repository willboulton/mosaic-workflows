#!/bin/bash

# $1 - base directory for sample, $2 - threads

# module load ruby/3.3.5
# module load R/4.2.1
# conda activate bio

SAMPLE=$(basename $1)

for f in $(ls $HOME/all_mags_322 | grep $SAMPLE); do 
    tar -xzOf $HOME/all_mags_322/$f ${f%.tar.gz}.fna > $1/bins/${f%.tar.gz}.fna
done

$HOME/software/DAS_Tool/src/Fasta_to_Contig2Bin.sh -e fna -i $1/bins > $1/das_c2b.tsv

seqkit grep -f <(cut -f 1 $1/das_c2b.tsv) $1/$SAMPLE.a.fna > $1/$SAMPLE.bins.fna

sed -e 's|_| |2' $1/$SAMPLE.a.faa |\
    seqkit grep -f <(cut -f 1 $1/das_c2b.tsv) |\
    awk 'BEGIN { OFS="_" }
    /^>/ {
        split($0, arr, " ")
        contig_id = arr[1]
        if (contig_id != prev_contig) {
            prev_contig = contig_id
            gene_number = 1
        } else {
            gene_number++
        }
        print contig_id "_" gene_number
        next
    } { print }' > $1/$SAMPLE.bins.faa

$HOME/software/DAS_Tool/DAS_Tool \
    -i $1/das_c2b.tsv \
    -c $1/$SAMPLE.bins.fna \
    -p $1/$SAMPLE.bins.faa \
    -o $1/bins_das/das \
    --threads $2 --score_threshold 0.0 --write_bins
