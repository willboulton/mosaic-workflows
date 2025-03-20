#!/bin/bash

# requires:
#  - singularity loaded 
#  - hmmer loaded
#  - python available

MC=/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue
EA=$MC/annotations/euk
f=$1

# merge braker with metaeuk using agat 
if [ -f br3/$f/braker.gff3 ]; then
    singularity exec ~/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif \
        agat_sp_complement_annotations.pl \
            -r $EA/rrna/$f/$f.gff \
            -a $EA/trna/$f/$f.gff \
            -a $EA/br3/$f/braker.gff3 \
            -a $EA/me/$f/$f.metaeuk.gff3 \
            --out $EA/pre/gff_struc_v1/$f.gff3
else:
    singularity exec ~/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif \
        agat_sp_complement_annotations.pl \
            -r $EA/rrna/$f/$f.gff \
            -a $EA/trna/$f/$f.gff \
            -a $EA/me/$f/$f.gff3 \
            --out $EA/pre/gff_struc_v1/$f.gff3
fi

# fix IDs in the gff3 - normalise final field
~/my_code/mosaic-workflows/workflow/scripts/euk_normalise_gff.py $EA $f

# get protein fasta files
singularity exec ~/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif \
    agat_sp_extract_sequences.pl \
    -g $EA/pre/gff_struc_v2/$f.gff3 \
    -f $MC/euk/$f.fna \
    -o $EA/faa/$f.faa --aa

# remove trailing stop characters as IPScan doesn't like it
sed -i -e 's|\*$||g' -e 's| .*||' $EA/faa/$f.faa

# run pfam scan
~/software/pfam_scan/pfam_scan.py -cpu 1 \
    -out $EA/pre/pfam/$f.pfam.csv \
    $EA/faa/$f.faa 
    ~/software/pfam_34_db

# merge pfam annotations into gff file 
# generate .pfam.txt file
#~/my_code/mosaic-workflows/arrays/pfam_add.py euk \
#    $EA/pre/pfam/$f.pfam.csv \
#    $EA/pre/gff_struc_v2/$f.gff3 \
#    $EA
