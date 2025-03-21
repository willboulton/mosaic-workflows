#!/usr/bin/env python3

from sys import argv, exit, stdout
from glob import glob
from functools import cache

from tqdm import tqdm 
import pandas as pd


PFAM_HEADER = [
    "gene_id",
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

PFAM_SCAN_TO_PFAM_HEADER_DICT = {
    "seq_id":"gene_id",
    "aln_start":"query_start", 
    "aln_end":"query_end", 
    "hmm_start":"subj_start", 
    "hmm_end":"subj_end", 
    "score":"bit_score"
}

PILOT_SAMPLES = {
    "3300049783",
    "3300046534",
    "3300047667",
    "3300047666",
    "3300046449",
    "3300045790",
    "3300045789",
    "3300046103",
    "3300046451",
    "3300047668",
    "3300046532",
    "3300046450",
    "3300047669",
    "3300045738",
    "3300046467",
}

@cache
def load_pfam_data_for_pftxt(pdat, add_col = None):
    
    pf = pd.read_csv(pdat, sep = ',')
    pf = pf.rename(columns = PFAM_SCAN_TO_PFAM_HEADER_DICT)
    pf['pfam_id'] = pf.hmm_acc.apply(lambda x: x.replace("PF", "pfam").split(".")[0])    
    pf['align_length'] = (pf.query_end - pf.query_start + 1).apply(int)
    pf['percent_identity'] = 100.0

    if add_col:
        pf['sample_id'] = str(add_col)

    return pf
   

def concat_pfam_data(ipf, pdat, opf, mode = 'sid_prefix'):

    pf = load_pfam_data_for_pftxt(pdat)

    if ipf == None:
        df = pf[PFAM_HEADER]
    else:
        df = pd.concat((
            pd.read_csv(
                ipf, 
                sep = '\t', 
                header = None, 
                names = PFAM_HEADER
            ), 
            pf[PFAM_HEADER]
        ))

    if mode == "ss_prefix':
        start_idx = 3
    elif mode == 'ss_noprefix':
        start_idx = 2
    elif mode == 'coa':
        start_idx = 3
    else:
        print("calls to this function must pass in an appropriate mode")
        exit(1)
    
    df['seq_id'] = df.gene_id.apply(lambda x: "_".join(x.split("_")[:start_idx]))
    df['gene_start'] = df.gene_id.apply(lambda x: int(x.split("_")[start_idx]))
    df['gene_end'] = df.gene_id.apply(lambda x: int(x.split("_")[start_idx + 1]))

    df = df.sort_values([
        "seq_id",
        "gene_start",
        "gene_end", 
        "query_start", 
        "query_end", 
        "evalue"
    ])

    df['bit_score'] = df['bit_score'].apply(str)
    df['percent_identity'] = df['percent_identity'].apply(str)

    df[PFAM_HEADER].to_csv(
        opf,
        sep = '\t', 
        index = False, 
        header = False, 
        float_format = '%.2e'
    ) 




def append_pfam_data_inplace(ipf, pf, mode):

    if mode == 'vir':

        df = pd.read_csv(
            ipf, 
            sep = '\t', 
            header = None, 
            names = PFAM_HEADER + ['contig_id']
        ) 

        df['sample_id'] = df.gene_id.apply(lambda x: str(x.split("_")[0]))
        df['gene_start'] = df.gene_id.apply(lambda x: int(x.split("_")[3]))
        df['gene_end'] = df.gene_id.apply(lambda x: int(x.split("_")[4]))

        if not set(df.sample_id).intersection(PILOT_SAMPLES):
            return

        if not set(df.gene_id).intersection(set(pf.gene_id)):
            return

        p2 = pf[pf.gene_id.isin(df.gene_id)]

    elif mode == "prok":
        
        df = pd.read_csv(
            ipf,
            sep = '\t',
            header = None,
            names = PFAM_HEADER
        )

        df['sample_id'] = ipf.split("/")[-1][:10]
        df['gene_start'] = df.gene_id.apply(lambda x: int(x.split("_")[2]))
        df['gene_end'] = df.gene_id.apply(lambda x: int(x.split("_")[3]))

        if not set(df.sample_id).intersection(PILOT_SAMPLES):
            return

        if not set(df.gene_id).intersection(set(pf.gene_id)):
            return

        p2 = pf[pf.gene_id.isin(df.gene_id) & pf.sample_id.isin(df.sample_id)]

    elif mode == "euk":
        exit(0)

    else:
        exit(1)
    
    df = pd.concat((df, p2))        

    df = df.sort_values([
        "contig_id",
        "gene_start",
        "gene_end", 
        "query_start", 
        "query_end", 
        "evalue"
    ])

    df['bit_score'] = df['bit_score'].apply(str)
    df['percent_identity'] = df['percent_identity'].apply(str)

    if mode == "vir":
        header = PFAM_HEADER + ['contig_id']
    if mode == "prok":
        header = PFAM_HEADER 
    if mode == "euk":
        header = PFAM_HEADER

    df[header].to_csv(
        ipf,
        sep = '\t', 
        index = False, 
        header = False, 
        float_format = '%.2e'
    ) 




def preprocess_gff(gff_file):

    df_gff = pd.read_csv(gff_file, sep = '\t', names = GFF_HEADER)
    df_gff['attributes_dict'] = df_gff.attributes.apply(
            lambda x: dict(tuple(y.split("=") for y in x.split(";")))
    )

    df_gff['gene_id'] = df_gff.attributes_dict.apply(lambda d: d.get("ID", "no_id"))

    return df_gff


def concat_pfam_to_gff(ipf, pdat, opf):
    
    df = preprocess_gff(ipf)
    pf = pd.read_csv(pdat, sep = ',').rename(columns = {"seq_id":"pfam_seq_id"})
    pf['pfam_id'] = pf.hmm_acc.apply(lambda x: x.split(".")[0])    

    df_edited = df[df.gene_id.isin(pf.pfam_seq_id)]
    #print(df_edited.shape)
    idx = df_edited.index
    
    df_edited = pd.merge(
        df_edited, 
        pf[['pfam_seq_id', 'pfam_id']], 
        left_on = 'gene_id', 
        right_on = 'pfam_seq_id', 
        how = 'inner'
    )

    #print(df_edited.shape)

    df_edited = df_edited.groupby('gene_id').agg({
        "seq_id":"first", 
        "source":"first",
        "type":"first", 
        "start":"first", 
        "end":"first", 
        "score":"first", 
        "strand":"first", 
        "phase":"first",
        "attributes":"first",
        "attributes_dict":"first", 
        "pfam_seq_id":"first", 
        "pfam_id":lambda x: ",".join(x)
   })

    #print(df_edited.shape) 
    #print(df_edited.head().T)

    def merge_pfam_into_dict(row):
        pf_result = row.attributes_dict.get('pfam', "")
        pf_result = set(pf_result.split(","))
        pf_result = pf_result.union(set(row.pfam_id.split(",")))
        pf_result = ",".join(pf_result).strip(",")
        
        row.attributes_dict['pfam'] = pf_result
        return row.attributes_dict


    df_edited['attributes_dict'] = df_edited.apply(merge_pfam_into_dict, axis = 1)

    df_edited['attributes'] = df_edited.attributes_dict.apply(
        lambda x: ";".join(f"{k}={v}" for k,v in x.items())
    )
    df_edited.index = idx

    #print(df_edited.index[:5])
    #print(df_edited.head().T)
    
    df = df[GFF_HEADER]
    df.update(df_edited[GFF_HEADER])

    df.to_csv(
        opf,
        sep = '\t',
        index = False,
        header = False
    )


if __name__ == "__main__":

    if argv[1] == "vir":

        ip_folder_path, pfam_folder_path = argv[2:]

        pfam_files = glob(f"{pfam_folder_path}/*")
        vir_files = glob(f"{ip_folder_path}/*")

        pf = pd.concat(
            (load_pfam_data_for_pftxt(f, add_col = f.split('/')[-1].split(".")[0]) for f in pfam_files), 
            axis = 0
        )
        pf['gene_id'] = pf.sample_id + "_" + pf.gene_id                
        pf['contig_id'] = pf.gene_id.apply(lambda x: "_".join(x.split("_")[:3]))
        pf['gene_start'] = pf.gene_id.apply(lambda x: int(x.split("_")[3]))
        pf['gene_end'] = pf.gene_id.apply(lambda x: int(x.split("_")[4]))

        #print(pf.head().T)

        if all(x.endswith(".pfam.txt") for x in vir_files):
            
            for f in tqdm(vir_files):
                append_pfam_data_inplace(f, pf, argv[1])       
            
            exit(0)

        elif all(x.endswith(".gff") for x in vir_files):
            exit(0)   

        else:
            pass

    elif argv[1] == "prok":
 
        ip_folder_path, pfam_folder_path = argv[2:]

        pfam_files = glob(f"{pfam_folder_path}/*")
        prok_files = glob(f"{ip_folder_path}/*")

        pf = pd.concat(
            (load_pfam_data_for_pftxt(f, add_col = f.split('/')[-1][:10]) for f in pfam_files), 
            axis = 0
        )
        pf['contig_id'] = pf.gene_id.apply(lambda x: "_".join(x.split("_")[:2]))
        pf['gene_start'] = pf.gene_id.apply(lambda x: int(x.split("_")[2]))
        pf['gene_end'] = pf.gene_id.apply(lambda x: int(x.split("_")[3]))

        #print(pf.head().T)

        if all(x.endswith(".pfam.txt") for x in prok_files):
            
            for f in tqdm(prok_files):

                print(f.split("/")[-1][:10])
                if f.split("/")[-1][:10] in PILOT_SAMPLES:
                    append_pfam_data_inplace(f, pf, argv[1])       
               
            exit(0)
 
        elif all(x.endswith(".gff") for x in prok_files):
            exit(0)   

        else:
            pass


    elif argv[1] == "euk":

        # euk mode is different - only takes one file at a time
        # the output is a fixed existing folder relative to the 
        # EA (euk annotation) folder given. 
        pfam_path, gff_path, ea_folder = argv[2:]

        mag_id = pfam_path.split("/")[-1].rstrip(".pfam.csv")

        if mag_id.find("33000") >= 0:
            mode = 'ss_prefix'
        else:
            mode = 'coa'

        concat_pfam_data(
            None,
            pfam_path, 
            f"{ea_folder}/pfam/{mag_id}.pfam.txt", 
            mode = mode
        )

        concat_pfam_to_gff(
            gff_path, 
            pfam_path, 
            f"{ea_folder}/gff/{mag_id}.gff"
        )
        
        exit(0)

    else:
        print("ERROR")
        print("Wasn't expecting this file extension - not sure what to do. ")
        exit(1)


