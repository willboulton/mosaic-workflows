#!/usr/env/bin python
import pandas as pd
import ete3
from sys import argv
from collections import Counter


def weighted_sum(counter_aggregation, node):
    
    # remove root node if it's not the only thing there
    # as this only drags estimates to be less specific
    
    if 1 in counter_aggregation:
        root_count = counter_aggregation.pop(1)
    
    if len(counter_aggregation) == 0:
        return root_count * node.get_distance(taxid_dict[1])
    
    return sum((w * node.get_distance(taxid_dict[n]) for n, w in counter_aggregation.items()))


def lineage(node):

    if node.up:
        x = lineage(node.up)
        x.append(node)
        return x

    return [node]

def minimise_weighted_sum(counter_aggregation):
    
    ans = {}
    
    def agg_weights(counter_aggregation, node):

        ans[node.name] = weighted_sum(counter_aggregation, node)

        if node.up:
            agg_weights(counter_aggregation, node.up)
    
    for n in counter_aggregation.keys():
        
        ans = {}
        agg_weights(counter_aggregation, taxid_dict[n])
    
    #print(ans)
    min_dist = min([d for n,d in ans.items()])
    #print("min dist is", min_dist)
    
    candidate_nodes = [taxid_dict[n] for n, d in ans.items() if d == min_dist]
    #print(candidate_nodes)
    
    if len(candidate_nodes) == 1:
        return candidate_nodes[0]
    
    lineages = set()
    for c in candidate_nodes:
        x = set(lineage(c))
        x.remove(c)
        lineages = lineages.union(x)

    #print("lineages:", lineages)

    cn = []
    for c in candidate_nodes:
        if not c in lineages:
            cn.append(c)

    #print(cn)

    if len(cn) == 1:
        return cn[0]

    return root.get_common_ancestor(cn)


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
        return "Viruses"
    return "other"


if __name__ == "__main__":

    # mmseqs input file
    mmseqs_tsv = argv[1]

    # STEP 1: First read in the NCBI taxonomy tree

    taxdump_folder = argv[2]
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


    # STEP 2: Next, read in mmseqs tsv file as snakemake input
    dft = pd.read_csv(
        input_mmseqs_tsv, 
        sep = '\t', 
        names = ['gene_id', 'taxid', 'tax_level', 'taxon', 'lineage']
    )
    
    dft['contig_id'] = dft.gene_id.apply(lambda x: "_".join(x.split("_")[:2]))
    

    dft = dft[['contig_id', 'taxid']]
    dft['taxid'] = dft.taxid.apply(int)
    dft = dft[dft.taxid > 1]
    dft = dft[dft.taxid.isin(taxid_dict)]
    
    dft = dft.groupby(
            'contig_id',
            sort=False,
            dropna=True).agg(
        lambda x: minimise_weighted_sum(Counter(x)).name
    )

    dft['contig_id'] = dft.index
    
    ans = pd.merge(
        dft, merged_names[['taxid', 'taxon', 'name']], on = 'taxid'
    )
    name_map = dict(zip(merged_names.taxid, merged_names.name))
    ans['new_lin'] = ans.taxid.apply(lambda x: ";".join([name_map.get(int(n.name), "") for n in lineage(taxid_dict[x])[2:]]))

    ans['superkingdom'] = ans.new_lin.apply(lambda x: fix_name(safe_get(0)(x)))

    ans = ans[['contig_id', 'taxid', 'superkingdom', 'taxon', 'name']]
    ans.to_csv(argv[3], sep = '\t', index = False)