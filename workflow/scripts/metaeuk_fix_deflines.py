#!/usr/bin/env python3 

import pandas as pd
import re
import itertools

from sys import argv, stdout
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# USAGE:

# Either ~/metaeuk_fix_deflines.py --inplace <basename> 
# Or ~/metaeuk_fix_deflines.py <basename> to write to stdout. 

# I'm expecting that you already ran metaeuk easy-annotate
# or whatever it's called, to get an .fas file. 
# Furthermore, you already ran metaeuk_to_gff3.py to generate
# a gff3 in the same folder. 

# This script will take the basename (without the extension)
# and generate a new <basename>.metaeuk.faa amino acid file
# from the <basename>.fas file given as input. 


BASE_NAME = argv[-1]

INFILE = f"{BASE_NAME}.fas"
GFF_FILE = f"{BASE_NAME}.gff3"

if argv[1] == "--inplace":
    OUTFILE = f"{BASE_NAME}.metaeuk.faa"


GFF_HEADER = [
    "contig_id",
    "software",
    "seqtype",
    "start",
    "end",
    "length",
    "strand",
    "dot",
    "desc"
]

# There are two kinds of header for the contigs, single sample, and coassembly
RE_PATTERN = r"\|\d+_Ga\d+_\d+"
RE_PATTERN += r"|\|[EIS]\d+_\d+_k\d+_\d+"


SEEN_DEFLINES={}

def parse_defline(defline):

    header = re.split(f"({RE_PATTERN})", defline)
    
    new_header = list(itertools.chain.from_iterable([x.strip("|").split("|") for x in header[1:]]))

    contig_id, gene_start, gene_end = new_header[0], new_header[5], new_header[6]

    result = f"{contig_id}_{int(gene_start) + 1}_{int(gene_end) + 1}"

    multiplicity = SEEN_DEFLINES.get(result, 0)

    if multiplicity == 0:
        
        SEEN_DEFLINES[result] = 1
        return f">{result}\n"

    else:
        multiplicity += 1
        SEEN_DEFLINES[result] = multiplicity
        return f">{result}_{multiplicity - 1}\n"


if __name__ == "__main__":

    result_str = ""
    
    with open(INFILE, 'r') as ipf:
        for line in ipf:
            if line[0] == ">":
                result_str += parse_defline(line)
            else:
                result_str += line

    if argv[1] == "--inplace":

        with open(OUTFILE, 'w') as opf:
            opf.write(result_str)

    else:
        stdout.write(result_str)
