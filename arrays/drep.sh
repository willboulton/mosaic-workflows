#!/bin/bash

#SBATCH --mem=256G             # memory
#SBATCH --ntasks=32           # n tasks=n threads
#SBATCH -p compute-64-512     # partition 
#SBATCH --time=00-23:59       # time (DD-HH:MM)
#SBATCH --job-name=drep_p       # Job name
#SBATCH -o drep_p.out       # Standard output log
#SBATCH -e drep_p.err       # Standard error log

#run the application

cd /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/

module load python/anaconda/2020.11/3.8
module load mash/2.2
module load mummer/4.0.0
module load fastANI/1.33
module load CheckM/1.1.2 

source /gpfs/home/qcf21tku/start_conda
conda activate bio

# alternative code was this:
#dRep dereplicate \
#    drep_bins_99 \
#    -p 16 \
#    -g drep_bins_input/*.fna \
#    --S_algorithm fastANI \
#    --S_ani 0.99 \
#    --ignoreGenomeQuality

dRep dereplicate \
    prok_drep_99 \
    -p 16 \
    -g prok/*.fna \
    --S_algorithm fastANI \
    --S_ani 0.99 \
    --completeness 50 \
    --contamination 10 \
    --genomeInfo prok_checkm/quals.csv
    