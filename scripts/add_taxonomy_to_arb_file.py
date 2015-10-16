#!/usr/bin/env python
import argparse


def add_tax_to_arbfile(taxonomy_file, arb_file):

    map_gtdb_id_to_tax = {}

    tax_fh = open(taxonomy_file, "rb")
    for line in tax_fh:
        splitline = line.rstrip().split("\t")
        map_gtdb_id_to_tax[splitline[0]] = splitline[1]
    tax_fh.close()

    current_id = None
    arb_fh = open(arb_file, "rb")
    for line in arb_fh:
        line = line.rstrip()
        splitline = line.split("=")
        if line[0:3] == "END":
            current_id = None
        elif splitline[0] == "db_name":
            current_id = splitline[1]
        elif (splitline[0] == "warning") and (current_id is not None):
            try:
                print "genome_tree_tax_string=%s" % map_gtdb_id_to_tax[current_id]
            except KeyError:
                pass
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
