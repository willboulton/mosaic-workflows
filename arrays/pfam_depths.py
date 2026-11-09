#!/usr/bin/env python3

import pandas as pd
import sys

PFAM_HEADER = [
 "seq_id", 
 "pfam_id",
 "percent_identity",
 "query_start",
 "query_end",
 "subj_start",
 "subj_end",
 "evalue",
 "bit_score",
 "align_length"
]

COG_HEADER = [

]

def main():

    # arguments: 
    # 1: the sampleid (33000*****)
    # 2: depth file
    # 3: pfam file

    with open(sys.argv[2], 'r') as depth_stream:
        df_depth = pd.read_csv(
            depth_stream, 
            sep = '\t'
        ).rename(
            columns = {"AvgFold":"depth", "ID":"contig_id"}
        )

  
    with open(sys.argv[3], 'r') as pfam_stream:
        df_pfam = pd.read_csv(
            pfam_stream, 
            sep = '\t', 
            header = None, 
            names = PFAM_HEADER
        )
        df_pfam['contig_id'] = df_pfam.seq_id.apply(lambda x: "_".join(x.split("_")[:2]))


    df_out = pd.merge(
                df_depth, 
                df_pfam, 
                on = 'contig_id', 
                how = 'inner'
    )[['pfam_id', 'depth']]

    tpm_scale = df_out.depth.sum() / (10 ** 6)


    df_out.groupby(
        'pfam_id'
    ).sum().divide(
        tpm_scale
    ).reset_index(
    ).to_csv(
        sys.stdout, 
        sep = '\t', 
        index = False, 
        header = ["pfam_id", str(sys.argv[1])]
    )


if __name__ == "__main__":

    main()
