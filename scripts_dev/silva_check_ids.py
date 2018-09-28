#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__prog_name__ = 'silva_check_ids.py'
__prog_desc__ = 'Check INSDC primary accession numbers for GTDB SSU sequences with those in SILVA.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2018'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import gzip
import logging
import ntpath
import argparse
from collections import defaultdict

import biolib.seq_io as seq_io


class Check(object):
    """Check INSDC primary accession numbers for GTDB SSU sequences with those in SILVA."""

    def __init__(self):
        """Initialize."""

        pass

    def run(self, ssu_gtdb_taxonomy_file, silva_parc_fasta_file):
        """Check INSDC primary accession numbers for GTDB SSU sequences with those in SILVA."""
        
        print('Reading SILVA INSDC accession numbers.')
        silva_ids = set()
        for seq_id, seq in seq_io.read_seq(silva_parc_fasta_file):
            silva_ids.add(seq_id)
        print('Read %d accession numbers.' % len(silva_ids))
            
        print('Checking GTDB SSU INSDC accession numbers.')
        missing_silva_acc = 0
        num_genes = 0
        for line in open(ssu_gtdb_taxonomy_file):
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            gene_id = line_split[1]
            gene_id = gene_id[0:gene_id.rfind('.')]
            start = int(line_split[2])
            stop = int(line_split[3])
            accession = '%s.%d.%d' % (gene_id, start, stop)
            num_genes += 1
            
            if accession not in silva_ids and (stop-start) > 300:
                print('Missing INSDC accession in SILVA for genome %s: %s (len=%d)' % (gid, accession, stop-start))
                missing_silva_acc += 1
                
        print('Identified %d of %d (%.2f%%) genes without a SILVA accession.' % (missing_silva_acc, num_genes, missing_silva_acc*100.0/num_genes))
 
if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ssu_gtdb_taxonomy_file', help='file with GTDB taxonomy for 16S rRNA genes specified by INSDC accession numbers')
    parser.add_argument('silva_parc_fasta_file', help="file with SILVA Parc sequences specified by INSDC accession numbers")

    args = parser.parse_args()

    try:
        p = Check()
        p.run(args.ssu_gtdb_taxonomy_file, args.silva_parc_fasta_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
