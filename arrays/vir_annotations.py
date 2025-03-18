#!/usr/bin/env python3

from tqdm import tqdm
import pandas as pd
import sys

BINNING_DIR = "/gpfs/data/mock_lab/will/mosaic_tmp/binning_dir"
VIRAL_OUT = "/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/annotations/vir"

GFF_DTYPES = {
    "seq_id": str, 
    "source": str,
    "type": str, 
    "start": int, 
    "end": int, 
    "score": str, 
    "strand": str, 
    "phase": str, 
    "attributes": str
}

GFF_HEADER = list(GFF_DTYPES.keys())


PFAM_DTYPES = {
    "gene_id":str,
    "pfam_id":str,
    "percent_identity":float,
    "query_start":int,
    "query_end":int,
    "subj_start":int,
    "subj_end":int,
    "evalue":float,
    "bit_score":float,
    "align_length":int
}

PFAM_HEADER = list(PFAM_DTYPES.keys())

def main(extn):
    
    # read in the viral bins
    df_vir = pd.read_csv(
        "/gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/mag_catalogue_to_bin.tsv", 
        sep = '\t', 
        header = None, 
        names = ['bin_id', 'contig_id']
    )
    df_vir = df_vir[df_vir.bin_id.str.startswith("vir")]
    df_vir['sample_id'] = df_vir.contig_id.apply(lambda x: x.split("_")[0])
    df_vir['contig_id_2'] = df_vir.contig_id.apply(lambda x: "_".join(x.split("_")[1:]))

    print(df_vir.head())


    def make_gff(sid):
        df_out = pd.read_csv(
            f"{BINNING_DIR}/{sid}/{sid}.a.gff", 
            sep = '\t', 
            header = None, 
            names = GFF_HEADER, 
            dtype = GFF_DTYPES, 
            comment = '#'
        )
        df_out = df_out[df_out.seq_id.isin(set(df_vir[df_vir.sample_id == sid].contig_id_2))]
        df_out['seq_id'] = df_out.seq_id.apply(lambda x: str(sid) + "_" + x)
        return df_out

    def make_pfam_tsv(sid):
        df_out = pd.read_csv(
            f"{BINNING_DIR}/{sid}/{sid}.a.pfam.txt", 
            sep = '\t', 
            header = None, 
            names = PFAM_HEADER, 
            dtype = PFAM_DTYPES, 
            comment = '#'
        )
        df_out['gene_id'] = df_out.gene_id.apply(lambda x: str(sid) + "_" + x)
        df_out['seq_id'] = df_out.gene_id.apply(lambda x: "_".join(x.split("_")[:3]))
        df_out = df_out[df_out.seq_id.isin(set(df_vir[df_vir.sample_id == sid].contig_id))]
        return df_out

    # make viral contig gff / tsv
    print("making viral gff/tsv file, this might take some time")
    if extn == "gff":
        extn2 = "gff"
        df = pd.concat(tqdm(make_gff(x) for x in set(df_vir.sample_id)))
    elif extn == "pfam":
        extn2 = "pfam.txt"
        df = pd.concat(tqdm(make_pfam_tsv(x) for x in set(df_vir.sample_id)))
    else:
        print("ERROR - must have one of gff or pfam in main")
        exit()
    print("finished viral gff/tsv")
    
    print(df.head())

    for f in tqdm(set(df_vir.bin_id)):
        contigs = set(df_vir[df_vir.bin_id == f].contig_id)
        df[df.seq_id.isin(contigs)].to_csv(
            f"{VIRAL_OUT}/{extn}/{f}.{extn2}", sep = '\t', index = False, header = None
        )


if __name__ == "__main__":

    main("pfam")


