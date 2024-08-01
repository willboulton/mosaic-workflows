#!/usr/bin/env python

from sys import argv, exit

contig_file = argv[1]
# assumes that the mapping lines are ordered the same
# as the contig ids in the contig file. 
mapping_file = argv[2]
output_file = argv[3]

with open(contig_file, 'r') as fh_in:
    with open(mapping_file, 'r') as fh_m:
        with open(output_file, 'w') as fh_out:
            
            scaffold, contig = tuple(next(fh_m).split("\t"))
            
            for l_in in fh_in:
                
                # id header lines in fasta file
                if l_in[0] == ">":
                    
                    newline = l_in.replace(contig, scaffold)
                    fh_out.writeline(newline)
                    
                    if newline == l_in:
                        exit("ERROR: Contig id did not match scaffold id - exiting")
                    
                    scaffold, contig = tuple(next(fh_m).split("\t"))

                # data line
                else:
                    fh_out.writeline(l_in)

print("Completed successfully")