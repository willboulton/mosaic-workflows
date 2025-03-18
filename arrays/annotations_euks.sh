#!/usr/bin/bash -e
#SBATCH --partition compute-64-512
#SBATCH --mem 32G
#SBATCH --time 720
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --array 1-1

module add RepeatMasker/4.1.2-p1
module add singularity/3.6.0
module add samtools/1.21
module add hmmer/3.3

source ~/start_conda
conda activate genemark

# Define paths
GENOME="your_genome.fasta"
OUTDIR="annotation_output"
METAOUT="metaeuk_out"
BRAKERDIR="braker_out"
PFAMDB="/path/to/Pfam-A.hmm"




# fix deflines from metaeuk

# run repeatmasker, softmasking
RepeatMasker -xsmall -species Algae -gff $GENOME

# if available, extract RNA-seq data from the contigs

# run braker, or Augustus

# run tRNA-scan
tRNAscan-SE -E --thread 8 --quiet \ 
    --gff ~/scratch/euglenozoa.1.split_1_v1_3300045738_braker3_test/trnascan_test.gff \
    ~/scratch/euglenozoa.1.split_1_v1_3300045738_braker3_test/genome.fa.masked

# run barrnap with euk settings
conda deactivate
conda activate bio


# make sure to generate an AA file with reasonable deflines

# run pfam_scan
~/software/pfam_scan/pfam_scan.py \
    $PFAMDB \
    
# incorperate pfam_scan results into a gff3
# Concat all annotations into a gff3 using AGAT

# tidy up intermediate files
