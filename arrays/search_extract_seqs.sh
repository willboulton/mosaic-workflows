###### Compress and concatenate files and generate databases ######

# Zip sequence files
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 24000 -t 1800 -n 1 --wrap "cd /gpfs/home/qcf21tku/scratch/Mosaic;gzip concat_contigs.fna && rm concat_contigs.fna"

# Concatenate fna files
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 24000 -t 240 -n 16 --wrap "cd /gpfs/data/mock_lab/will/mosaic_tmp/img_files_12_3_2024;ls | cut -f 1 -d '.' | ~/software/parallel/parallel-20240222/bin/parallel -j 16 \"zcat {}.tar.gz | tar xfO - {}/{}.a.fna | sed 's/>/>{}_/g'\" > ~/scratch/Mosaic/concat_contigs.fna"

# Concatenate faa files
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 24000 -t 240 -n 16 --wrap "cd /gpfs/data/mock_lab/will/mosaic_tmp/img_files_12_3_2024;ls | cut -f 1 -d '.' | ~/software/parallel/parallel-20240222/bin/parallel -j 16 \"zcat {}.tar.gz | tar xfO - {}/{}.a.faa | sed 's/>/>{}_/g'\" > ~/scratch/Mosaic/concat_contigs.faa"

# Generate database files
sbatch -p hmem-4T -c 1 --mem 800000 -t 600 -n 32 --qos hmem --wrap "cd /gpfs/home/qcf21tku/scratch/Mosaic;module load ncbi-blast/2.10.1+;makeblastdb -in concat_contigs.fna -dbtype nucl -parse_seqids -out mosaic_contigs.fna.blastdb -title mosaic_contigs.fna.blastdb"

sbatch -p hmem-4T -c 1 --mem 800000 -t 600 -n 32 --qos hmem --wrap "cd /gpfs/home/qcf21tku/scratch/Mosaic;module load ncbi-blast/2.10.1+;makeblastdb -in concat_contigs.faa -dbtype prot -parse_seqids -out mosaic_contigs.faa.blastdb -title mosaic_contigs.faa.blastdb"


###### Sequence searches - HMMER, table lookup from PFam (parallel) and BLAST ######

# Run HMMs 
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 16000 -t 360 -n 1 --wrap "module add hmmer/3.3;cd ~/hmms;hmmsearch -E 1e-5 --tblout ~/scratch/blast_results/big2_vs_mosaic_contigs_faa_1em5.hmmtblout round2_big2.hmm ~/scratch/Mosaic/concat_contigs.faa.gz"

# Get all sequence ids based on PFAM domains
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 16000 -t 240 -n 24 --wrap "cd /gpfs/data/mock_lab/will/mosaic_tmp/img_files_12_3_2024;ls | cut -f 1 -d '.' | ~/software/parallel/parallel-20240222/bin/parallel -j 24 \"zcat {}.tar.gz | tar xfO - {}/{}.a.pfam.txt | grep -wFf ~/query_seqs/crispr_pfams.txt\" > ~/scratch/blast_results/crispr_pfam_seqs.tsv"

# Blast - blastp
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 2000 -t 60 -n 8 --wrap "module add ncbi-blast/2.10.1+;blastp -query ~/query_seqs/aureochrome.faa -db ~/scratch/Mosaic/mosaic_contigs_faa_blastdb/mosaic_contigs.faa.blastdb -out ~/scratch/blast_results/aureochrome_vs_mosaic_faa.v3.blast.out -seg no -evalue 1e-5 -num_threads 8 -outfmt \"6 qseqid sseqid pident length qlen slen qcovs evalue bitscore\""

# Blast - blastn
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 4000 -t 120 -n 8 --wrap "module add ncbi-blast/2.10.1+;blastn -query ~/query_seqs/big1.fna -db ~/scratch/Mosaic/mosaic_contigs_fna_blastdb/mosaic_contigs.fna.blastdb -out ~/scratch/blast_results/big1_vs_mosaic_faa.v1.blast.out -dust no -evalue 1e-5 -num_threads 8 -outfmt \"6 qseqid sseqid pident length qlen slen qcovs evalue bitscore\""

###### Extract sequences from concatenated faa and fna files ######

# Extract sequences from file
sbatch -p compute-64-512,compute-24-96 -c 1 --mem 8000 -t 120 -n 8 --wrap "source ~/start_conda;conda activate blobology;mytmp=\$(mktemp);cat ~/scratch/blast_results/big2_vs_mosaic_contigs_faa_1em5.hmmtblout | grep -P \"^33000\d+\_Ga\d+\_\d+\_\d+\_\d+\" -o > \${mytmp};seqkit grep -j 8 -o ~/scratch/blast_results/big2_vs_mosaic_contigs_faa_1em5.hmmresults.faa -f \$mytmp /gpfs/data/mock_lab/will/concat_contigs.faa.gz;rm -f \$mytmp"

sbatch -p compute-64-512,compute-24-96 -c 1 --mem 8000 -t 120 -n 8 --wrap "source ~/start_conda;conda activate blobology;mytmp=\$(mktemp);cat ~/scratch/blast_results/aureochrome_vs_mosaic_faa.v2.blast.out | cut -f 2 > \$mytmp;seqkit grep -j 8 -o ~/scratch/blast_results/aureochrome_vs_mosaic_contigs_faa_1em5.blast.out.faa -f \$mytmp /gpfs/data/mock_lab/will/concat_contigs.faa.gz;rm -f \$mytmp"
