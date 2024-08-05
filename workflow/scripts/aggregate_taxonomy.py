#!/usr/bin/env python 

from sys import argv
import pandas as pd

def coarse_classify(sk):
    
    if sk in ["Eukaryota", "Virus", "Viruses", "eukarya"]:
        return "EUV"
    if sk in ['Unknown', 'unknown', 'organelle']:
        return "UNK"
    if sk in ['Archaea', "Bacteria", "prokarya", "archaea", "bacteria"]:
        return "PRO"
    


def combine_estimates(row):

    # Centrifuge is the worst so in the event of any ties, 
    # we drop that one first. 
    # Hence, we can first check if any of the other 2 are equal 
    # (and not Unknown). 
    
    # If the two best methods agree
    if (row.gp_phylo == row.gp_mmseqs):
        if row.gp_phylo != "UNK":
            return row.gp_phylo
        else:
            if row.gp_tiara != "UNK":
                return row.gp_tiara
            else:
                return row.gp_cent
    
    # Otherwise phylo != mmseqs, use tiara as a tiebreak
    else: 
        if (row.gp_phylo == row.gp_tiara):
            if row.gp_phylo == "UNK":
                return row.gp_mmseqs
            else:
                return row.gp_phylo
        
        elif (row.gp_mmseqs == row.gp_tiara):
            if row.gp_mmseqs == "UNK":
                return row.gp_phylo
            else:
                return row.gp_mmseqs   

        # phylo, mmseqs, tiara all different
        else: 
            return row.gp_cent
            

if __name__ == "__main__":

    # read phylodist
    print("reading phylodist")
    df_phylo = pd.read_csv(argv[1], 
        sep = '\t', 
        names = ["gene_id", 'hit1', 'hit2', 'pc_id', "lineage"]
    )
    
    df_phylo['contig_id'] = df_phylo.gene_id.apply(lambda x: "_".join(x.split("_")[:2]))
    df_phylo['superkingdom_phylo'] = df_phylo.lineage.apply(lambda x: x.split(";")[0])
    df_phylo = df_phylo[['contig_id', 'superkingdom_phylo']].groupby('contig_id').agg(
        lambda x: pd.Series.mode(x)[0]
    ).reset_index()

    print(df_phylo.head())
    print(len(df_phylo))

    # read mmseqs
    print("reading mmseqs")
    df_mmseqs = pd.read_csv(
        argv[2],
        names = ['contig_id', 'taxid', 'superkingdom_mmseqs', 'lineage', 'taxon', 'name'], 
        header = 0, 
        sep = '\t'
    )
    df_mmseqs = df_mmseqs[['contig_id', 'superkingdom_mmseqs']]

    print(df_mmseqs.head())
    print(len(df_mmseqs))
    
    # read tiara
    print("reading tiara")
    df_tiara = pd.read_csv(
        argv[3], 
        sep = '\t'
    ).rename(columns = {
        'sequence_id':'contig_id', 
        'class_fst_stage': 'superkingdom_tiara',
        'class_snd_stage': 'dk'
    })[['contig_id', 'superkingdom_tiara']]

    print(df_tiara.head())
    print(len(df_tiara))
  

    # read centrifuge
    print("reading centrifuge")
    df_cent = pd.read_csv(
        argv[4], 
        sep = '\t'
    )

    taxa_superkingdom = pd.read_csv(f"{argv[5]}/ncbi_superkingdom_2024-03-10.csv")

    taxa_superkingdom = taxa_superkingdom.fillna("Unknown")

    print("combining centrifuge with taxonomy")
    df_cent = pd.merge(
        df_cent,
        taxa_superkingdom,
        left_on = 'taxID', 
        right_on = 'tax_id', 
        how = 'left'
    ).fillna("Unknown")

    
    df_cent = df_cent[['readID', 'superkingdom']].rename(
        columns = {'readID':'contig_id', 'superkingdom':'superkingdom_cent'}
    )

    print(df_cent.head())
    print(len(df_cent))

    print("merging dataframes")
    df_total = pd.merge(df_phylo, df_mmseqs, on = 'contig_id', how = 'outer').fillna("Unknown")
    df_total = pd.merge(df_total, df_tiara, on = 'contig_id', how = 'outer').fillna("Unknown")
    df_total = pd.merge(df_total, df_cent, on = 'contig_id', how = 'outer').fillna("Unknown")

    df_total['gp_phylo'] = df_total.superkingdom_phylo.apply(coarse_classify)
    df_total['gp_mmseqs'] = df_total.superkingdom_mmseqs.apply(coarse_classify)
    df_total['gp_tiara'] = df_total.superkingdom_tiara.apply(coarse_classify)
    df_total['gp_cent'] = df_total.superkingdom_cent.apply(coarse_classify)


    print("getting concensus kingdom from all tools")
    df_total['gp_concensus'] = df_total.apply(combine_estimates, axis = 1)

    print(df_total.head())
    print(len(df_total))
    # write EVU contigs to file
    print("saving output")
    df_total[df_total.gp_concensus == "EVU"][['contig_id']].to_csv(argv[6], header = None)