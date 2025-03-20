#!/bin/bash
#SBATCH --mem=64G                    # memory
#SBATCH --ntasks=8
#SBATCH -p compute-64-512,compute-24-96            # partition
#SBATCH --time=00-23:59             # time (DD-HH:MM)
#SBATCH --job-name=merge_genes     # Job name
#SBATCH --array=1-13
#SBATCH -o merge_genes-%j.out       # Standard output log
#SBATCH -e merge_genes-%j.err       # Standard error log

#run the application
module add singularity/3.6.0

conda activate bio

MC=/mag_catalogue
EA=$MC/annotations/euk
f=$1

# merge braker with metaeuk using agat 
if [ -f br3/$f/braker.gff3 ]; then 
    singularity exec ~/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif \
        agat_sp_complement_annotations.pl \
            -r $EA/rrna/$f/$f.gff \
            -a $EA/trna/$f/$f.gff \
            -a $EA/br3/$f/braker.gff3 \
            -a $EA/me/$f/$f.gff3 \
            --out $EA/pre/gff_struct_v1/$f.gff3
else:
    singularity exec ~/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif \
        agat_sp_complement_annotations.pl \
            -r $EA/rrna/$f/$f.gff \
            -a $EA/trna/$f/$f.gff \
            -a $EA/me/$f/$f.gff3 \
            --out $EA/pre/gff_v1_struct/$f.gff3
fi

# fix IDs in the gff3 - normalise final field

# get protein fasta files
singularity exec ~/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif \
    agat_sp_extract_sequences.pl \
    -g $EA/pre/gff_v1_struct/$f.gff3 \
    -f $MC/euk/$f.fna \
    -o $EA/faa/$f.faa --aa

# fix deflines in faa files 
# remove trailing stop characters as IPScan doesn't like it

# run interproscan or pfam_scan

# merge pfam_scan results into gff3