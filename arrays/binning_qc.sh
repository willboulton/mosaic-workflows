# Eukcc
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 32000 -t 180 -n 8 --wrap "module add singularity/3.6.0;singularity run ~/software/eukcc/eukcc_latest.sif folder --db ~/software/eukcc/eukcc2_db_ver_1.2 --out comebin_eukcc --suffix .fa --threads 8 comebin_out"

# Checkm2
swrap -p compute-64-512 --mem 64G -c 1 -n 24 -t 600 --wrap "source ~/start_conda; conda activate checkm2; checkm2 predict --threads 24 --input /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/prok --output-directory /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/prok_checkm --force -x fna --database_path /gpfs/home/qcf21tku/scratch/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"

# Move bins
cd /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir
for g in $(ls | grep 33000); do 
    for f in $(awk '$10 > 49 && $11 < 10' ${g}/bins_das/das_DASTool_summary.tsv | cut -f 1); do 
        seqkit grep -f <(grep -w $f $g/bins_das/das_DASTool_contig2bin.tsv | cut -f 1) $g/${g}.a.fna > ${g}/${g}_${f}.fa
    done
done


# Generate graph file for euk deduplications
for f in `cat euk_samples.txt`; do 
    for a in `ls euk_with_duplicates/ | grep $f`; do 
        for b in `ls euk_with_duplicates/ | grep $f`; do 
            printf "${a}\t${b}\t"
            comm -12 \
                <(grep '>' euk_with_duplicates/$a | tr -d '>' | sort) \
                <(grep '>' euk_with_duplicates/$b | tr -d '>' | sort) | wc -l
        done
    done
done > euk_deduplications.txt