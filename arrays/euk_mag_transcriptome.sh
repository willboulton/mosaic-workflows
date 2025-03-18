#!/bin/bash -e
#SBATCH --partition compute-64-512,compute-24-96
#SBATCH --mem 24G
#SBATCH -t 180
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --array 3-356

module add samtools/1.16.1

MC=/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue
EUKT=$MC/annotations/euk/transcriptome

a=( $(head -n $SLURM_ARRAY_TASK_ID /gpfs/home/qcf21tku/euks_all.txt | tail -n 1) )
f=${a[0]}

for i in $(seq 1 14); do 
    mkdir -p $EUKT/$f/m$i
    N=$(( $i * 10 ))
    for bam in $(ls $MC/transcriptome_mappings/*.bam | head -n $N | tail -n 10); do 
        len=$(samtools view -h $bam | head -n 2 | tail -n 1 | cut -f 3 -d ":")
        samtools view --regions-file \
            <(grep ">" $MC/euk/$f.fna | tr -d '>' | sed "s|$|\t0\t$len|g") \
            -o $EUKT/$f/m$i/${bam##*/} $bam
    done
    samtools merge -o $EUKT/$f/m$i/merged.bam $EUKT/$f/m$i/*.bam
done

mkdir -p $EUKT/$f/m15
for bam in $(ls $MC/transcriptome_mappings/*.bam | tail -n 7); do 
    len=$(samtools view -h $bam | head -n 2 | tail -n 1 | cut -f 3 -d ":")
    samtools view --regions-file \
        <(grep ">" $MC/euk/$f.fna | tr -d '>' | sed "s|$|\t0\t$len|g") \
        -o $EUKT/$f/m15/${bam##*/} $bam
done
samtools merge -o $EUKT/$f/m15/merged.bam $EUKT/$f/m15/*.bam

find $EUKT/$f/m*/*.bam -empty -delete

mkdir -p $EUKT/$f/all
samtools merge -o $EUKT/$f/all/all_merged.bam $EUKT/$f/m*/merged.bam && \
    rm -rf $EUKT/$f/m*


