#!/bin/bash
#SBATCH --mem=128G                    # memory
#SBATCH -p compute-64-512           # partition 
#SBATCH --time=00-23:59             # time (DD-HH:MM)
#SBATCH --job-name=will_gtdb_tk     # Job name
#SBATCH -o gtdb_tk_job-%j.out       # Standard output log
#SBATCH -e gtdb_tk_job-%j.err       # Standard error log
#SBATCH --ntasks=24                 # n tasks=n threads

module load Prodigal/2.6.3
module load hmmer/3.3
module load mash/2.2
module load FastTree/2.1.11
module load pplacer/1.1.alpha19


source /gpfs/home/qcf21tku/start_conda
conda activate gtdb_tk

export GTDBTK_DATA_PATH=/gpfs/home/qcf21tku/scratch/gtdb_release220

BASEPATH=/gpfs/scratch/qcf21tku/Mosaic/jgi_havoc_pilot_mags_zipped

gtdbtk classify_wf \
    --genome_dir $BASEPATH/jgi_mags_genomes \
    --out_dir $BASEPATH/gtdb_tk_out \
    --mash_db $BASEPATH/mash/mash.msh \
    --cpus 24 \
    --pplacer_cpus 24