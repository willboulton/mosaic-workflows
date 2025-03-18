#!/bin/bash

mkdir -p $1/all_bins_multi_auto

for g in `ls $1/metabat_multi`; do
    for f in `grep -Po ">33000\d\d\d\d\d" $1/metabat_multi/$g | sort | uniq | tr -d '>'`; do 
        seqkit grep -rp $f $1/metabat_multi/$g > $1/all_bins_multi_auto/${g%.fa}_$f.fna
    done
done


for g in `ls $1/vamb_multi/bins`; do
    for f in `grep -Po ">33000\d\d\d\d\d" $1/vamb_multi/bins/$g | sort | uniq | tr -d '>'`; do 
        seqkit grep -rp $f $1/vamb_multi/bins/$g > $1/all_bins_multi_auto/${g%.fna}_$f.fna
    done
done

