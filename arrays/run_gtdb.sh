#!/bin/bash
#SBATCH --mem=256G                    # memory
#SBATCH -p compute-64-512           # partition 
#SBATCH --time=01-23:59             # time (DD-HH:MM)
#SBATCH --job-name=will_gtdb_tk     # Job name
#SBATCH -o gtdb_tk_job-%j.out       # Standard output log
#SBATCH -e gtdb_tk_job-%j.err       # Standard error log
#SBATCH --ntasks=48                 # n tasks=n threads

module load Prodigal/2.6.3
module load hmmer/3.3
module load mash/2.2
module load FastTree/2.1.11
module load pplacer/1.1.alpha19


source /gpfs/home/qcf21tku/start_conda
conda activate gtdb_tk

export GTDBTK_DATA_PATH=/gpfs/home/qcf21tku/scratch/gtdb_release220

BASEPATH=/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue

gtdbtk classify_wf \
    --genome_dir $BASEPATH/prok \
    --out_dir $BASEPATH/prok_gtdb_tk_out \
    --mash_db $BASEPATH/prok_mash/mash.msh \
    --cpus 48 \
    --pplacer_cpus 48
