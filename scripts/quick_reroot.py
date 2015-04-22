#!/usr/bin/env python
import argparse
import os
import sys
import tempfile
import subprocess

import ete2

def reroot(tree_filepath, output_prefix=None):
    t = ete2.Tree(tree_filepath)

    ancestor = t.get_common_ancestor("IMG_2264867063", "IMG_2524023065")
    t.set_outgroup(ancestor)

    if output_prefix is not None:
       t.write(format=0, outfile=output_prefix)
    else:
       print t.write(format=0)
 
if __name__ == '__main__':

    # create the top-level parser
    parser = argparse.ArgumentParser(prog='quick_reroot.py')
    parser.add_argument('-t', dest='tree', required=True,
                        help='Newick tree to reroot.')
    parser.add_argument('-o', dest='output_prefix', required=True, 
                        help='Output prefix for tax to tree')
    args = parser.parse_args()

    sys.stderr.write("\nThis script uses ETE. If you use it, consider citing the following:\nJaime Huerta-Cepas, Joaquin Dopazo and Toni Gabaldon.\nETE: a python Environment for Tree Exploration. BMC Bioinformatics 2010, 11:24.\n\n")
    sys.stderr.flush()
    
    reroot(args.tree, args.output_prefix) 
