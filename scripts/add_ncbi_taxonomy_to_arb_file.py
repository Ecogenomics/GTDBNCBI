#!/usr/bin/env python
import argparse
import os
import sys


def add_tax_to_arbfile(taxonomy_file, arb_file):
    
    map_organism_name_to_tax = {}
    
    tax_fh = open(taxonomy_file, "rb")
    for line in tax_fh:
        splitline = line.rstrip().split("\t")
        organism_name = splitline[1]
        tax_string = splitline[0]
        split_organism_name = organism_name.split()
        if (len(split_organism_name) < 2):
            continue
        
        organism_name = " ".join(split_organism_name[0:2])
        
        if organism_name in map_organism_name_to_tax:
            old_tax_string = map_organism_name_to_tax[organism_name]
            if old_tax_string != tax_string:
                if old_tax_string is not None:
                    sys.stderr.write(
                        "%s has multiple tax strings. Will not be able to assign\n    Old:%s\n    New:%s\n" %
                        (organism_name, old_tax_string, tax_string)
                    )
                    sys.stderr.flush()
                map_organism_name_to_tax[organism_name] = None
            continue
        map_organism_name_to_tax[organism_name] = tax_string
    tax_fh.close()
    
    current_organism = None
    current_tax = None
    arb_fh = open(arb_file, "rb")
    for line in arb_fh:
        line = line.rstrip()
        splitline = line.split("=")
        #sys.stderr.write(str(splitline) + "\n")
        if line[0:3] == "END":
            current_organism = None
        
        elif splitline[0] == "organism":
            current_organism = splitline[1]
            matches = set()
            for tax_assoc_organism_name in map_organism_name_to_tax:
                if current_organism.startswith(tax_assoc_organism_name):
                    tax_to_add = map_organism_name_to_tax[tax_assoc_organism_name]
                    if tax_to_add is not None:
                        matches.add(tax_to_add)
                    #sys.stderr.write("%s: Adding %s - %s\n" % (current_organism, tax_assoc_organism_name, map_organism_name_to_tax[tax_assoc_organism_name]))
            if len(matches) > 1:
                sys.stderr.write("Too many tax hits for %s:\n    %s\n" % (current_organism, "\n    ".join(matches)))
                current_organism = None
            elif len(matches) == 1:
                #sys.stderr.write("MATCH!!! %s for %s\n" % (str(matches), current_organism))
                current_tax = matches.pop()
            else:
                current_organism = None
                
        elif (splitline[0] == "warning") and (current_organism is not None):
            print "ncbi_tax_string=%s" % current_tax
        print line
    arb_fh.close()
 
if __name__ == '__main__':

    # create the top-level parser
    parser = argparse.ArgumentParser(prog='quick_reroot.py')
    parser.add_argument('-t', dest='taxfile', required=True,
                        help='Newick tree to reroot.')
    parser.add_argument('-a', dest='arbfile', required=True, 
                        help='Output prefix for tax to tree')
    args = parser.parse_args()
    
    add_tax_to_arbfile(args.taxfile, args.arbfile) 