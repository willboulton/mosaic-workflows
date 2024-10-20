#!/usr/bin/env python

from sys import exit
from scipy import sparse
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import umap
import subprocess
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

    parser.add_argument('--fasta_in', type = str, default = '', 
                        help = "Input fasta file of all contigs - only necessary if you're using fasta_out."
                       )
    
    parser.add_argument('--fasta_out', type = str, default = '', 
                        help = "Generate fasta files of the clusters (bins) in this directory." + 
                               "If left empty, doesn't produce the bins."
                       )

    parser.add_argument('--contig_map', type = str, default = '', 
                        help = 'Mapping file in case you want to rename scaffold ids to contigs ids.')

    parser.add_argument('--contignames', type = str, default = '', 
                        help = 'File containing the contig ids in the same order as they appear in the embeddings file.')

    parser.add_argument('--umap_only', action = 'store_true', help = 'Only do the umap bit.')

    
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
        print("Loading npz matrix")
        mat = np.load(args.embeddings).get("arr_0")
        print("VAMB Latent space loaded")
    elif args.embeddings_format == "COMEBin":
        print("Loading COMEBin latent space")
        df_emb = pd.read_csv(args.embeddings, sep = '\t', index_col = 0)
        mat = df_emb.to_numpy()
        contignames = pd.DataFrame(df_emb.index, columns = ['contig_id'])
        print("COMEBin Latent space loaded")
    else:
        exit("unknown file, only VAMB npz, or COMEBin tsv allowed")

    print("running umap - default parameters")
    data = umap.UMAP().fit(mat)

    print("Saving transformed output coordinates file")
    coords = data.transform(mat)
    extn = args.umap_out.split(".")[-1]
    if extn == "npy":
        np.save(
            args.umap_out[:-4], 
            coords
        )

    if args.embeddings_format == "VAMB" and not args.contignames:
        print("VAMB embedding matrix doesn't provide contig names, you need --contignames.")
        exit()
    
    # provided by vamb
    if args.contignames:
        contignames = pd.read_csv(args.contignames, names = ['contig_id'])
    
    if extn in ["tsv", "txt"]:
        df_coords = pd.DataFrame(coords, columns = ['umap_0', 'umap_1'])
        df_coords['contig_id'] = contignames.contig_id
        df_coords = df_coords[['contig_id', 'umap_0', 'umap_1']]
        df_coords.to_csv(args.umap_out, sep = '\t', index = False)
        print(f"wrote {df_coords.shape} dataframe to {args.umap_out}")
    elif extn == "csv":
        df_coords = pd.DataFrame(coords, columns = ['umap_0', 'umap_1'])
        df_coords['contig_id'] = contignames.contig_id
        df_coords = df_coords[['contig_id', 'umap_0', 'umap_1']]
        df_coords.to_csv(args.umap_out, sep = ',', index = False)
    else:
        exit(f"The parameter --umap_out, {args.umap_out} must have an extension one of npy, tsv, txt, or csv")

    if args.umap_only:
        exit()
        
    # get the umap weighted graph and put it in a readable form
    # (to, from) value
    print("generating ijv graph")
    umap_graph = data.graph_.tocoo()
    graph_coo = sparse.find(umap_graph)
    
    # put the graph matrix into ijv form
    ijv_mat = np.vstack(graph_coo).T
    
    # Optionally, transform the 3rd column between 0 to 10
    transform_between_0_10 = args.transform
    
    if transform_between_0_10:
        ijv_mat[:, 2] = 10 + np.log10(ijv_mat[:, 2] + 0.0000000001)
    
    # the 'value' is the umap likelihood of connection. 
    print("saving ijv graph to temporary file")
    with NamedTemporaryFile() as tmpf1:
        with NamedTemporaryFile() as tmpf2:
            np.savetxt(
                f"{tmpf1.name}", 
                ijv_mat, 
                delimiter = '\t', 
                fmt = ['%i', '%i', '%.10e']
            )

            # run MCL algorithm like this:

            print("running mcl: " + f"mcl {args.umap_out}.ijv_graph.tsv --abc -I {args.inflation} -o {tmpf2.name}")
            subprocess.run([
                "mcl", 
                f"{tmpf1.name}", 
                "--abc", 
                "-I", f"{args.inflation}", 
                "-o", f"{tmpf2.name}"
            ])

            # -I inflation parameter is most important, smaller gives larger cluster but is slower
            # Try experimenting with mclcm heirarchical clutering

            
            with open(tmpf2.name, 'r') as mcl_results:
                mcl_clusters = [
                    [int(x) for x in cluster.split("\t") if x.strip() != ""] for cluster in mcl_results.read().split("\n")
                ]
            
            #contig_ids = pd.read_csv(args.contig_map, header = False).to_list()
            mcl_clusters = [
                [contignames.contig_id.iloc[x] for x in cluster] for cluster in mcl_clusters if len(cluster) > 0
            ]

            with open(f"{args.mcl_out}", 'w') as mcl_clusters_out:
                
                mcl_clusters_out.write(f"mcl_bin_id\tcontig_id\n")
                
                for idx, cluster in enumerate(mcl_clusters):                
                    for c in cluster:
                        mcl_clusters_out.write(f"mcl.{idx + 1}\t{c}\n")
                    
                    if args.fasta_out and args.fasta_in:

                        with open(f"{args.fasta_out}/mcl.{idx + 1}.fna", 'w') as seq_fna:
                            with NamedTemporaryFile(mode = "w") as tmpf:
                                
                                tmpf.write("\n".join(cluster))
                                subprocess.run([
                                    "seqkit", "grep", "--quiet", "-nf", tmpf.name, args.fasta_in
                                ], stdout = seq_fna)

