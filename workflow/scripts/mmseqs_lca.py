#!/usr/bin/env python

import pandas as pd
import ete3
from sys import argv
from functools import lru_cache
from collections import Counter

@lru_cache(maxsize = None)
def get_lineage(node):

    if node.up:
        x = get_lineage(node.up)
        x.append(node)
        return x

    return [node]


def safe_get(n, sep=';'):
    
    def func(string):
        x = string.split(sep)
        if len(x) <= n:
            return ""
        return x[n]
    
    return func

def fix_name(name_str):
    if name_str in ["Archaea", "Bacteria", "Eukaryota"]:
        return name_str
    if "vir" in name_str:
        return "Virus"
    return "Unknown"


def generate_lineage_string(taxid, name_map):

    lineage = get_lineage(taxid_dict[taxid])
    
    start = min(len(lineage) - 1, 2)
    
    list_of_taxa = [name_map.get(int(n.name), "") for n in lineage[start:]]
    
    return ";".join(list_of_taxa)


def ncbi_taxonomy_data(taxdump_folder):
    df = pd.read_csv(
        f"{taxdump_folder}/nodes.dmp", 
        sep = "|", 
        header = None
    )[[0,1,2]].tail(-1)
    
    df.columns = ['taxid', 'parent_taxid', 'taxon']
    
    names = pd.read_csv(
        f"{taxdump_folder}/names.dmp", 
        sep = '|',
        names = ['taxid', 'name', 'other_name', 'name_type', 'something']
    ).applymap(lambda x: x.strip() if type(x) == str else x)
    
    names = names[names.name_type == "scientific name"]
    
    merged_names = pd.merge(df, names, on = 'taxid')
    merged_names['taxon'] = merged_names.taxon.apply(lambda x: x.strip())
    
    root = ete3.Tree(name = 1)
    taxid_dict = { 1: root }
    
    for t in df.itertuples():
        taxid, parentid = t[1], t[2]
        
        if not taxid_dict.get(parentid):
            taxid_dict[parentid] = ete3.Tree(name = parentid)
        
        if not taxid_dict.get(taxid):
            taxid_dict[taxid] = ete3.Tree(name = taxid)
            
        taxid_dict[parentid].add_child(taxid_dict[taxid])

    return merged_names, taxid_dict

class TaxData():
    __slots__ = ("leaf", "weight", "depth")
    
    def __init__(self, leaf = False, weight = 0, depth = 0):
        
        self.leaf = leaf
        self.weight = weight
        self.depth = depth
        
    def __repr__(self):
        return f"({self.leaf}, {self.weight}, {self.depth})"

    
uninformative_taxa = [
    0, # unclassified
    1, # root node
    131567, # cellular organisms
    2787823, # unknown 1
    2787854 # unknown 2
]

superkingdoms = [
    2, # bacteria
    10239, # viruses
    2157, # archaea
    2759 # eukaryota
]

def aggregate_tax_estimates(counter_aggregation):

    global done
    done += 1
    if (done % 10000 == 0):
        print(f"\t{done} of {n_annotated_contigs}")

    if len(counter_aggregation) == 1:
        return list(counter_aggregation.keys())[0]

    total_weight = sum(w for idx, w in counter_aggregation.items() if not (idx in superkingdoms or not idx in taxid_dict))
    ignored_weight = sum(w for idx, w in counter_aggregation.items() if not ((not idx in superkingdoms) or not idx in taxid_dict))
    
    # if the ignored weight is too much, we include it again
    if ignored_weight > 5 * total_weight:
        total_weight += ignored_weight

    ans = {}

    for tax, weight in counter_aggregation.items():
        depth = len(get_lineage(taxid_dict.get(tax))) if tax in taxid_dict else 0
        ans[tax] = TaxData(leaf = True, weight = weight, depth = depth)

    for tax, weight in counter_aggregation.items():

        lineage_tax = [t.name for t in get_lineage(taxid_dict.get(tax))][:-1] if tax in taxid_dict else [1]

        for depth, l in enumerate(lineage_tax):
            if l in ans:
                ans[l].leaf = False
                ans[l].weight += weight
                ans[l].depth = max(ans[l].depth, depth)
            else:
                ans[l] = TaxData(leaf = False, weight = weight, depth = depth)

    
    leaves = { k:v for k,v in ans.items() if v.leaf }
    if len(leaves) == 0:
        print(leaves, ans)
        for a in ans:
            print([t.name for t in get_lineage(taxid_dict.get(a))])

    while True:

        maxdepth = max(t.depth for _,t in leaves.items())
        for l in leaves:
            if 2 * leaves[l].weight > total_weight:
                return l

        best_weight_so_far = -1
        best_tax_id = 1
        for l in leaves:
            if leaves[l].depth == maxdepth and leaves[l].weight > best_weight_so_far:
                best_tax_id = l

        test_tax = leaves.pop(best_tax_id)
        parent_id = get_lineage(taxid_dict[best_tax_id])[-2].name
        leaves[parent_id] = ans[parent_id]



if __name__ == "__main__":

    # mmseqs input file
    mmseqs_tsv = argv[1]

    # STEP 1: First read in the NCBI taxonomy tree
    taxdump_folder = argv[2]
    merged_names, taxid_dict = ncbi_taxonomy_data(taxdump_folder)
    print(f"read ncbi taxonomy: {len(merged_names)} rows")

    # STEP 2: Next, read in mmseqs tsv file as snakemake input
    dft = pd.read_csv(
        mmseqs_tsv, 
        sep = '\t', 
        names = ['gene_id', 'taxid', 'tax_level', 'taxon', 'lineage']
    )

    total_rows = len(dft)
    print(f"read mmseqs tsv file: {total_rows} rows")
    
    dft['contig_id'] = dft.gene_id.apply(lambda x: "_".join(x.split("_")[:2]))
    

    dft = dft[['contig_id', 'taxid']]
    dft['taxid'] = dft.taxid.apply(int)
    dft = dft[~dft.taxid.isin(uninformative_taxa)]
    dft = dft[dft.taxid.isin(taxid_dict)]
    filtered_rows = len(dft)
    
    print(f"removed taxids: {filtered_rows} rows")
    dft = dft.groupby(
        'contig_id',
        sort=False,
        dropna=True
    ).agg(Counter).reset_index()
    
    print("finished counter aggregations")
    n_annotated_contigs = len(dft)
    print(f"{n_annotated_contigs} contigs")

    done = 0
    dft['taxid'] = dft.taxid.apply(aggregate_tax_estimates)
    print("finished taxonomy aggregations")
    
    ans = pd.merge(
        dft, merged_names[['taxid', 'taxon', 'name']], on = 'taxid'
    )
    print("merged answer with ncbi taxonomy")
    name_map = dict(zip(merged_names.taxid, merged_names.name))
    ans['lineage'] = ans.taxid.apply(lambda x: generate_lineage_string(x, name_map))

    ans['superkingdom'] = ans.lineage.apply(lambda x: fix_name(safe_get(0)(x)))

    print("saving output file")
    ans = ans[['contig_id', 'taxid', 'superkingdom', 'lineage', 'taxon', 'name']]
    ans.to_csv(argv[3], sep = '\t', index = False)