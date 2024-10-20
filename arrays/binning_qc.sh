# Eukcc
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 32000 -t 180 -n 8 --wrap "module add singularity/3.6.0;singularity run ~/software/eukcc/eukcc_latest.sif folder --db ~/software/eukcc/eukcc2_db_ver_1.2 --out comebin_eukcc --suffix .fa --threads 8 comebin_out"

# Checkm2

