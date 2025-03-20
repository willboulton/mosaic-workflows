#!/usr/bin/env python3

from sys import argv
import pandas as pd

GFF_HEADER = [
    "seq_id", 
    "source",
    "type", 
    "start", 
    "end", 
    "score", 
    "strand", 
    "phase", 
    "attributes"
]

def preprocess_gff(gff_file):

    df_gff = pd.read_csv(gff_file, sep = '\t', names = GFF_HEADER, comment = '#')
    df_gff['attributes_dict'] = df_gff.attributes.apply(
            lambda x: dict(tuple(y.split("=") for y in x.split(";")))
    )

    return df_gff

ID_MAP = {}
ID_COUNTER = {}

def make_id(row):

    attr_dict = row.attributes_dict

    if row.source.startswith("barrnap"):
        return attr_dict
    if row.source.startswith("tRNAscan"):
        return attr_dict

    if "Parent" in attr_dict:
        attr_dict["Parent"] = ID_MAP[attr_dict["Parent"]]

    if row.source.startswith("metaEuk"):
        new_id = f"{row.seq_id}_{row.start}_{row.end}_{row.type}"

        if new_id in ID_COUNTER:
            count = ID_COUNTER.get(new_id)
            ID_COUNTER[new_id] = count + 1
            new_id = new_id + str(count)
        else:
            ID_COUNTER[new_id] = 1

    else:
        new_id f"{row.seq_id}_{row.start}_{row.end}_{attr_dict["ID"]}"
    
    ID_MAP[attr_dict["ID"]] = new_id
    attr_dict["ID"] = new_id
    return attr_dict

if __name__ == "__main__":

    euk_ann_basedir = argv[1].rstrip("/")
    sid = argv[2]

    gff_in = f"{euk_ann_basedir}/pre/gff_struc_v1/{sid}.gff3"
    gff_out = f"{euk_ann_basedir}/pre/gff_struc_v2/{sid}.gff3"

    df = preprocess_gff(gff_in)

    df['attributes_dict'] = df.apply(make_id, axis = 1)

    df['attributes'] = df.attributes_dict.apply(
        lambda x: ";".join(f"{k}={v}" for k,v in x.items())
    )
    
    df[GFF_HEADER].to_csv(
        gff_out,
        sep = '\t',
        header = None,
        index = None
    )