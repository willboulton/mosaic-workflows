#!/usr/env/python

import sys
import pandas as pd
import ete3
import csv
from functools import lru_cache
from mmseqs_lca import get_lineage, ncbi_taxonomy_data


@lru_cache(maxsize = 64000)
def classify_tax_id(taxid, taxid_to_parent):
    eukaryota_taxid = 2759
    virus_taxid = 10239
    bacteria_taxid = 2
    archaea_taxid = 2157
    
    lineage = set(x.name for x in get_lineage(taxid, taxid_to_parent))

    if bacteria_taxid in lineage:
        return 'Bacteria'
    elif eukaryota_taxid in lineage:
        return 'Eukaryota'
    elif virus_taxid in lineage:
        return 'Virus'
    elif archaea_taxid in lineage:
        return 'Archaea'
    else:
        return 'Unknown'

def filter_reads(centrifuge_output, taxonomy_dir, output_file):
    
    merged_names, taxid_dict = ncbi_taxonomy_data(taxdump_dir)
    
    with open(centrifuge_output, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        for row in reader:
            taxid = int(row['taxID'])
            classification = classify_tax_id(taxid, taxid_to_parent)
            if classification in {'Eukaryota', 'Virus', 'Unknown'}:
                writer.writerow([row['seqID']])

if __name__ == "__main__":
    centrifuge_output = sys.argv[1]
    taxonomy_dir = sys.argv[2]
    output_file = sys.argv[3]
    filter_reads(centrifuge_output, taxonomy_dir, output_file)
