#!/usr/bin/bash -e
#SBATCH --partition compute-64-512,compute-24-96
#SBATCH --mem 32G
#SBATCH --time 600
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --array 101-354

module add singularity/3.6.0

a=( $(head -n $SLURM_ARRAY_TASK_ID /gpfs/home/qcf21tku/euks_all.txt | tail -n 1) )

EA=/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/annotations/euk

GN=${a[0]}

[[ ! -f $EA/br3/$GN/braker.gff3 ]] && singularity exec -B ${PWD}:${PWD} \
    ~/software/braker3/braker3.sif braker.pl \
    --genome=$EA/repeatmasker/$GN/$GN.fna.masked \
    --min_contig=500 \
    --threads 8 \
    --gff3 --workingdir=$EA/br3/$GN/ \
    --esmode --species sp_unk_$SLURM_ARRAY_TASK_ID --overwrite

