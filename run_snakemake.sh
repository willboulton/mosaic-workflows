#!/usr/bin/env bash

#SBATCH --mem=16G                           # memory
#SBATCH --ntasks=1                          # n tasks=n threads
#SBATCH -p compute-64-512,compute-24-64     # partition 
#SBATCH --time=10-00:00                     # time (DD-HH:MM)
#SBATCH --job-name=snakemake-main           # Job name
#SBATCH -o snakemake-%j.out                 # Standard output log
#SBATCH -e snakemake-%j.err                 # Standard error log


. ~/.bashrc

cd /gpfs/home/qcf21tku/my_code/mosaic-workflows
conda activate snakemake
snakemake --profile profile_slurm/