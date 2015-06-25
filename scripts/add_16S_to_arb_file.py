#!/usr/bin/env python
import argparse
import os
import sys


def add_16S_blast_to_arbfile(blast_file, arb_file):
    
    map_gtdb_id_to_hits = {}
    
    blast_fh = open(blast_file, "rb")
    for line in blast_fh:
        splitline = line.rstrip().split("\t")
        (gtdb_id, accession, contig, hits, length, tax) = splitline
        if gtdb_id not in map_gtdb_id_to_hits:
            map_gtdb_id_to_hits[gtdb_id] = {}
        score = 0
        if tax == "-":
            tax = ""
        if accession in map_gtdb_id_to_hits[gtdb_id]:
            #if map_gtdb_id_to_hits[gtdb_id][accession][2] < score:
            map_gtdb_id_to_hits[gtdb_id][accession] = (hits, length, score, tax)
        else:
            map_gtdb_id_to_hits[gtdb_id][accession] = (hits, length, score, tax)
        
    blast_fh.close()
    
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
                print "blast_hits_16S_count=%i" % len(map_gtdb_id_to_hits[current_id].keys())
                print "blast_hits_16S=%s" % " # ".join(([ "%s %s/%s (%i%%) %s" % (accession, hits, length, (int(hits) * 100)/int(length), tax) for accession, (hits, length, score, tax) in map_gtdb_id_to_hits[current_id].items()]))
            except KeyError:
                pass
        print line
    arb_fh.close()
 
if __name__ == '__main__':

    # create the top-level parser
    parser = argparse.ArgumentParser(prog='add_16S_to_arb_file.py')
    parser.add_argument('-b', dest='blastfile', required=True,
                        help='Blast file.')
    parser.add_argument('-a', dest='arbfile', required=True, 
                        help='Arb file')
    args = parser.parse_args()
    
    add_16S_blast_to_arbfile(args.blastfile, args.arbfile) 