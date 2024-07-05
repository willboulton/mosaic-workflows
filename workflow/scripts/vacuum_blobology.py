#!/usr/env/python

from sys import argv, exit
import numpy as np
import pandas as pd
import umap
from scipy import sparse
import os
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description='Apply UMAP to VAMB or COMEBin embeddings, and run MCL clustering.')
    
    parser.add_argument('--embeddings', type=str, help='File to the embeddings.')

    parser.add_argument('--inflation', '-i', type=float, default = 1.2, help = 'MCL inflation parameter.')
    
    parser.add_argument('--umap_out', type=str, help = 'Path to UMAP 2d embeddings output file.')

    parser.add_argument('--mcl_out', type=str, help = 'Path to MCL cluster output file.')
    
    parser.add_argument('--embeddings_format', default='VAMB', choices = ["VAMB", "COMEBin"],
                        help='VAMB|COMEBin, depending on the format of the embeddings file (npz or tsv).')

    parser.add_argument('--transform', type = bool, default = False, 
                        help = 'log-transform similarities to between 0 and 10.')

    parser.add_argument('--fasta_out', type = str, default = '', 
                        help = "Generate fasta files of the clusters (bins) in this directory." + 
                               "If left empty, doesn't produce the bins."
                       )

    parser.add_argument('--contig_map', type = str, default = '', 
                        help = 'Mapping file in case you want to rename scaffold ids to contigs ids.')

    
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = parse_args()
    # $1 - latent embedding space - npz or tsv file
    # $2 - output file of umap embeddings - vamb_2k_latent_umap_default.npy
    # $3 - output file of ijv graph format for mcl clustering - vamb_umap_default_ijv_graph.tsv
    # $4 - MCL inflation parameter e.g. 1.2
    
    # load the latent space matrix (N x 32 size)
    
    if args.embeddings_format == "VAMB":
        mat = np.load(args.embeddings).get("arr_0")
    elif args.embeddings_format == "COMEBin":
        pass
    else:
        exit("unknown file extension, only npz, tsv, csv, or txt allowed")
    
    data = umap.UMAP().fit(mat)
    
    coords = data.transform(mat)
    np.save(
        args.umap_out, 
        coords
    )
    
    # get the umap weighted graph and put it in a readable form
    # (to, from) value
    umap_graph = data.graph_.tocoo()
    graph_coo = sparse.find(umap_graph)
    
    # put the graph matrix into ijv form
    ijv_mat = np.vstack(graph_coo).T
    
    # Optionally, transform the 3rd column between 0 to 10
    transform_between_0_10 = args.transform
    
    if transform_between_0_10:
        ijv_mat[:, 2] = 10 + np.log10(ijv_mat[:, 2] + 0.0000000001)
    
    # the 'value' is the umap likelihood of connection. 
    np.savetxt(
        f"{args.umap_out}.ijv_graph.tsv", 
        ijv_mat, 
        delimiter = '\t', 
        fmt = ['%i', '%i', '%.10e']
    )
    
    # run MCL algorithm like this:
    # mcl vamb_umap_default_ijv_graph_log_sim.tsv --abc -I 1.2
    
    os.subprocess(f"mcl {args.umap_out}.ijv_graph.tsv --abc -I {args.inflation}")
    # -I inflation parameter is most important, smaller gives larger cluster but is slower
    # Try experimenting with mclcm heirarchical clutering
    
    for 
    os.subprocess(f"seqkit grep -f {} > {argv[4]}")
