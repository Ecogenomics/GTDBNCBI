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

__prog_name__ = 'ncbi_extract_nt_genomes.py'
__prog_desc__ = 'Extract genes in nucleotide space from GenBank file.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2016'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import tempfile
from collections import defaultdict

from numpy import mean, std

class ExtractGenes(object):
  """Extract genes in nucleotide space."""

  def __init__(self):
    pass

  def run(self, genome_dir_file):

    for line in open(genome_dir_file):
        line_split = line.strip().split('\t')
        
        genome_id = line_split[0]
        genome_path = line_split[1]
        genome_dir_id = os.path.basename(os.path.normpath(genome_path))
        
        genbank_file = os.path.join(line_split[1], genome_dir_id + '_genomic.gbff')
        nt_gene_file = os.path.join(line_split[1], genome_dir_id + '_protein.fna')
        
        if not os.path.exists(nt_gene_file):
            cmd = './genbank2nt_fasta.pl %s > %s' % (genbank_file, nt_gene_file)
            print cmd
            os.system(cmd)
        else:
            print 'Nucleotide gene file already exists: %s' % genome_id

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_dir_file', help='file indicating path to all genome directories')
  
    args = parser.parse_args()

    try:
        p = ExtractGenes()
        p.run(args.genome_dir_file)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
