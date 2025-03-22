#!/usr/bin/env python3

from glob import glob 
import pandas as pd

PFAM_HEADER = [
    "seq_id", 
    "pfam_id", 
    "pc_iden", 
    "aln_start", 
    "aln_end", 
    "hmm_start", 
    "hmm_end", 
    "evalue", 
    "bit_score", 
    "aln_length"
]

pfam_files = glob("*/pfam/*")

df_pfams = pd.concat((
    pd.read_csv(
        f, 
        sep = '\t', 
        header = None, 
        names = PFAM_HEADER
    ).pfam_id.value_counts(
    ).to_frame(
    ).rename(
        columns = {
            "pfam_id" : f.split('/')[-1].strip(".gff"),
        }
    )
    for f in pfam_files), 
    axis = 1
).fillna(0).applymap(int).sort_index(
).T.to_csv(
    "/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/annotations/bin_pfam_counts.tsv",
    sep = '\t', 
    index = True, 
    header = True
)

