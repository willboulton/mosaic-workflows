#!/bin/bash -e
#SBATCH --partition compute-64-512,compute-24-96
#SBATCH --mem 80G
#SBATCH -t 360
#SBATCH -c 16
#SBATCH -n 1
#SBATCH --array 1-16

N=$(( $SLURM_ARRAY_TASK_ID * 2 ))

source ~/start_conda
conda activate bio
module load metaeuk
BASE=/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue

for q in $(head -n $N ~/euk_ss_remaining.txt | tail -n 2); do

    metaeuk easy-predict \
        $BASE/euk/$q.fna \
        /gpfs/data/mock_lab/reference-data/MOSAIC/MOSAIC_PROT_REF_DB/MOSAIC_NR_u50_o100_DB \
        $BASE/annotations/euk/me/$q \
        ~/scratch/${q}_metaeuk_tmp \
        --threads 16 --split-memory-limit 70G

    ~/my_code/mosaic-workflows/workflow/scripts/metaeuk_to_gff3.py \
        -f $BASE/annotations/euk/me/$q.fas \
        -g $BASE/annotations/euk/me/$q.gff3

    rm -rf ~/scratch/${q}_metaeuk_tmp

done
