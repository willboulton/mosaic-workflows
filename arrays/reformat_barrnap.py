#!/usr/bin/env python3

from sys import argv, stdin, stdout, exit
import pandas as pd

gff_names = [
    "contig_id", 
    "software", 
    "ann_type", 
    "start", 
    "end", 
    "evalue", 
    "strand",
    "dot", 
    "comment" 
]

def fix_contig_id(contig_id):

    if contig_id.startswith("33000"):
        contig_id = "_".join(contig_id.split("_")[1:])        
        
    return contig_id


def add_id_field(row):
   
    comment = row.comment
    comment = f"ID={row.contig_id}_{row.start}_{row.end};" + comment
    return comment

try: 
    df_barrnap = pd.read_csv(stdin, sep = '\t', comment = "#", header = None, names = gff_names)
    df_barrnap['contig_id'] = df_barrnap['contig_id'].apply(fix_contig_id)
    df_barrnap['comment'] = df_barrnap.apply(add_id_field, axis = 1)
    df_barrnap.to_csv(stdout, sep = '\t', index = False, header = False)

except ValueError:
    exit(0)

