#!/usr/bin/env bash

SAMPLE_ID=$(cut -f 5 /gpfs/data/mock_lab/will/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

tar -xzf {input} $SAMPLE_ID/$SAMPLE_ID.a.fna -O > {output.fna}
tar -xzf {input} $SAMPLE_ID/$SAMPLE_ID.a.faa -O > {output.faa}
tar -xzf {input} $SAMPLE_ID/$SAMPLE_ID.a.phylodist.txt -O > {output.phylodist}
tar -xzf {input} $SAMPLE_ID/$SAMPLE_ID.a.map.txt -O > {output.map}
tar -xzf {input} $SAMPLE_ID/$SAMPLE_ID.a.depth.txt -O > {output.depth}



module add centrifuge/1.0.4

centrifuge -x /gpfs/data/mock_lab/will/EUKFINDER_CENTRIFUGE/Centrifuge_NewDB_Sept2020 \
    -U /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/$SAMPLE_ID/$SAMPLE_ID.a.fna \
    -f -p {threads} -S {output}


source activate tiara
tiara -m 2000 -t 2 \
    -i /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/$SAMPLE_ID/$SAMPLE_ID.a.fna \
    -o /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID/$SAMPLE_ID.a.tiara.txt


/gpfs/home/qcf21tku/my_code/mosaic-workflows/workflow/scripts/run_mmseqs.sh \
    $SAMPLE_ID \
    /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID \
    /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID/$SAMPLE_ID_mmseqs \
    /gpfs/data/mock_lab/reference-data/MOSAIC/MOSAIC_PROT_REF_DB/MOSAIC_NR_u50_o100_DB \
    2 \
    /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/$SAMPLE_ID


/gpfs/home/qcf21tku/my_code/mosaic-workflows/workflow/scripts/mmseqs_lca.py \
    /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID/$SAMPLE_ID.faa.mms2.tsv \
    /gpfs/home/qcf21tku/software/taxdump /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID/$SAMPLE_ID.faa.mms2lca



/gpfs/home/qcf21tku/my_code/mosaic-workflows/workflow/scripts/aggregate_taxonomy.py \
    /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID/$SAMPLE_ID.a.phylodist.txt \
    /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID/$SAMPLE_ID.faa.mms2lca \
    /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID/$SAMPLE_ID.a.tiara.txt \
    /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID/$SAMPLE_ID.a.cf \
    /gpfs/home/qcf21tku/software/taxdump \
    /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/$SAMPLE_ID/$SAMPLE_ID.a.map.txt \
    /gpfs/data/mock_lab/will/mosaic_tmp/taxo_tmpdir/$SAMPLE_ID/$SAMPLE_ID.a.euks.txt
