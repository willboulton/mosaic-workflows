#!/bin/bash -e
#SBATCH --partition compute-64-512
#SBATCH --mem 256G
#SBATCH -t 360
#SBATCH -c 16
#SBATCH -n 1
#SBATCH --array 1-149

MC=/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue


a=( $(head -n $SLURM_ARRAY_TASK_ID /gpfs/home/qcf21tku/mosaic_metaT.txt | tail -n 1) )

f=${a[0]}
n=$(echo $f | cut -d/ -f 6)

source ~/start_conda
conda activate anvio
module add samtools/1.21

if [ ! -s $MC/transcriptome_mappings/$n.bam ]; then
    strobealign -t 16 $MC/mag_catalogue.fna \
        -r 150 --use-index --interleaved ${a[0]} > $MC/transcriptome_mappings/$n.bam
fi

if [ ! -s $MC/transcriptome_mappings/$n.sorted.bam ]; then
    samtools sort -m 4G -@ 16 -T $HOME \
        $MC/transcriptome_mappings/$n.bam > $MC/transcriptome_mappings/$n.sorted.bam && \
        rm $MC/transcriptome_mappings/$n.bam
fi