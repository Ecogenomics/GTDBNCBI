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

__prog_name__ = 'ltp_ssu_classification_fetch.py'
__prog_desc__ = 'Fetch BLAST results against the Living Tree Project (LTP)'

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
import ntpath
import argparse
from collections import defaultdict

from biolib.parallel import Parallel
from biolib.external.execute import check_dependencies


class LTP(object):
    """Fetch BLAST results against the Living Tree Project (LTP)."""

    def __init__(self):
        pass

    def run(self, ncbi_refseq_dir, ncbi_genbank_dir, user_genome_dir, output_file):
        """Fetch BLAST results against the Living Tree Project (LTP)."""

        fout = open(output_file, 'w')
        fout.write('accession\tquery_id\ttaxonomy\tlength\tsubject_id\tevalue\tbitscore\talign_len\tperc_identity\n')

        # generate metadata for NCBI assemblies
        for ncbi_genome_dir in [ncbi_refseq_dir, ncbi_genbank_dir]:
            print('Reading NCBI assembly directories in %s.' % ncbi_genome_dir)
            for domain in ['archaea', 'bacteria']:
                domain_dir = os.path.join(ncbi_genome_dir, domain)
                for species_dir in os.listdir(domain_dir):
                    full_species_dir = os.path.join(domain_dir, species_dir)
                    for assembly_dir in os.listdir(full_species_dir):
                        accession = assembly_dir[0:assembly_dir.find('_', 4)]
                        
                        ltp_blast_file = os.path.join(full_species_dir, assembly_dir, 'rna_ltp_132', 'ssu.taxonomy.tsv')
                        if os.path.exists(ltp_blast_file):
                            with open(ltp_blast_file) as f:
                                f.readline()
                                for line in f:
                                    fout.write('%s\t%s' % (accession, line))

        # generate metadata for user genomes
        if user_genome_dir != 'NONE':
            print('Reading user genome directories.')
            for user_id in os.listdir(user_genome_dir):
                full_user_dir = os.path.join(user_genome_dir, user_id)
                if not os.path.isdir(full_user_dir):
                    continue

                for genome_id in os.listdir(full_user_dir):
                    ltp_blast_file = os.path.join(full_user_dir, genome_id, 'rna_ltp_132', 'ssu.taxonomy.tsv')
                    if os.path.exists(ltp_blast_file):
                        with open(ltp_blast_file) as f:
                            f.readline()
                            for line in f:
                                fout.write('%s\t%s' % (accession, line))
                                
        fout.close()

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ncbi_refseq_dir', help='base directory leading to NCBI RefSeq archaeal and bacterial genome assemblies')
    parser.add_argument('ncbi_genbank_dir', help='base directory leading to NCBI RefSeq archaeal and bacterial genome assemblies')
    parser.add_argument('user_genome_dir', help='base directory leading to user genomes or NONE to skip')
    parser.add_argument('output_file', help='output file')

    args = parser.parse_args()

    try:
        p = LTP()
        p.run(args.ncbi_refseq_dir, args.ncbi_genbank_dir, args.user_genome_dir, args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
