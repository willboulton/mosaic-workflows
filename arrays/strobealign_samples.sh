sbatch -p compute-64-512 -c 1 --mem 32000 -t 360 -c 8 --array 1-1 --wrap "module add samtools/1.16.1;source ~/start_conda;conda activate anvio;SAMPLE=\$(ls /gpfs/data/mock_lab/will/fastq/tmp_fastq/*.fastq.gz | cut -f 8 -d '/' | cut -f 1 -d '-' | head -n \${SLURM_ARRAY_TASK_ID} | tail -1);mkdir -p /gpfs/home/qcf21tku/scratch/\${SAMPLE}/bam;strobealign -t 8 /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/\${SAMPLE}/\${SAMPLE}.a.fna --interleaved /gpfs/data/mock_lab/will/fastq/tmp_fastq/\${SAMPLE}-filter-METAG.fastq.gz | samtools sort -@ 8 -o /gpfs/home/qcf21tku/scratch/\${SAMPLE}/bam/\${SAMPLE}.sorted.bam"

swrap --partition compute-64-512 --mem 92G -t 180 -c 8 -n 1 --array 1-1 --array-file ~/mosaic_metaT.txt --dry-run --wrap "f=\${a[0]}; n=\$(echo \$f | cut -d/ -f 6); source ~/start_conda; conda activate anvio; module add samtools/1.21; strobealign -t 8 /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/mag_catalogue.fna --use-index --interleaved \${a[0]} | samtools sort -@ 8 -m 4G -T \${HOME} > /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/transcriptome_mappings/\$n.bam"


```
#!/bin/bash -e
#SBATCH --partition compute-64-512
#SBATCH --mem 256G
#SBATCH -c 16
#SBATCH -t 120

source ~/start_conda
conda activate anvio
cd /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/
strobealign --create-index mag_catalogue.fna -r 150
```

```
#!/bin/bash -e
#SBATCH --partition compute-64-512
#SBATCH --mem 256G
#SBATCH -t 60
#SBATCH -c 16
#SBATCH -n 1
#SBATCH --array 2-147

MC=/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue


a=( $(head -n $SLURM_ARRAY_TASK_ID /gpfs/home/qcf21tku/mosaic_metaT.txt | tail -n 1) )

f=${a[0]}
n=$(echo $f | cut -d/ -f 6)

source ~/start_conda
conda activate anvio
module add samtools/1.16.1


if [ ! -s $MC/transcriptome_mappings/$n.bam ]; then
    strobealign -t 16 $MC/mag_catalogue.fna \
        --use-index --interleaved ${a[0]} > $MC/transcriptome_mappings/$n.bam
fi

if [ ! -s $MC/transcriptome_mappings/$n.sorted.bam ]; then
    samtools sort -m 4G -@ 16 -T $HOME \
        $MC/transcriptome_mappings/$n.bam > $MC/transcriptome_mappings/$n.sorted.bam && \
        rm $MC/transcriptome_mappings/$n.bam
fi
```