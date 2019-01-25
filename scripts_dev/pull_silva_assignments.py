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

__prog_name__ = 'pull_silva_assignments.py'
__prog_desc__ = 'Create table with 16S rRNA assignments for GTDB genomes.'

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
import argparse


class PullHits(object):
    """Create table with 16S rRNA assignments for GTDB genomes."""

    def __init__(self):
        """Initialize."""
        
        self.ssu_dir = 'rna_silva'
        
    def run(self, 
                genome_file, 
                min_iden, 
                min_len,
                min_aln_len,
                output_file):
        """Create table with 16S rRNA assignments for GTDB genomes."""
        
        fout = open(output_file, 'w')
        
        bHeader = True
        for idx, line in enumerate(open(genome_file)):
            if idx % 1000 == 0:
                print('Processed %d genomes.' % idx)
                
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            gpath = line_split[1]
            
            rna_path = os.path.join(gpath, self.ssu_dir)
            taxonomy_file = os.path.join(rna_path, 'ssu.taxonomy.tsv')
            
            if not os.path.exists(taxonomy_file):
                continue
                
            with open(taxonomy_file) as f:
                headers = f.readline().strip().split('\t')
                if bHeader:
                    fout.write('Genome ID\t' + '\t'.join(headers) + '\n')
                    bHeader = False
                
                pident_index = headers.index('blast_perc_identity')
                aln_len_index = headers.index('blast_align_len')
                ssu_len_index = headers.index('length')
                
                for line in f:
                    line_split = line.strip().split('\t')
                    
                    pident = float(line_split[pident_index])
                    aln_len = int(line_split[aln_len_index])
                    ssu_len = int(line_split[ssu_len_index])
                    
                    if pident >= min_iden and aln_len >= min_aln_len and ssu_len >= min_len:
                        fout.write('%s\t%s' % (gid, line))
                        
        fout.close()

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_file', help='file indicating path to genome files')
    parser.add_argument('--min_iden', help='minimum identity to report assignment', default=0.95)
    parser.add_argument('--min_len', help='minimum length of 16S rRNA gene to report assignment', default=500)
    parser.add_argument('--min_aln_len', help='minimum alignment length to report assignment', default=500)
    parser.add_argument('output_file', help="output file")

    args = parser.parse_args()

    try:
        p = PullHits()
        p.run(args.genome_file, 
                args.min_iden,
                args.min_len,
                args.min_aln_len,
                args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
