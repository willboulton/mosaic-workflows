#!/bin/bash

#SBATCH --mem=256G             # memory
#SBATCH --ntasks=32           # n tasks=n threads
#SBATCH -p hmem-4T,hmem-754,compute-64-512     # partition
#SBATCH --qos hmem 
#SBATCH --time=01-23:59       # time (DD-HH:MM)
#SBATCH --job-name=drep_e3       # Job name
#SBATCH -o drep_e3.out       # Standard output log
#SBATCH -e drep_e3.err       # Standard error log

#run the application

cd /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/

module load python/anaconda/2020.11/3.8
module load mash/2.2
module load mummer/4.0.0
module load fastANI/1.33
module load CheckM/1.1.2 

source /gpfs/home/qcf21tku/start_conda
conda activate bio


# This was in there before - not needed if re-running (apparently)
#-g euk/*.fna \

dRep dereplicate \
    euk_drep_99 \
    -g euk/*.fna \
    -p 32 \
    --S_algorithm fastANI \
    --S_ani 0.99 \
    --completeness 30 \
    --contamination 15 \
    --genomeInfo euk_eukcc/euk_eukcc_quals.csv
