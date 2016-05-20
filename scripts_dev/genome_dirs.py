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

__prog_name__ = 'ncbi_genome_dirs.py'
__prog_desc__ = 'Produce file indicating the directory of each genome.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
from collections import defaultdict

class GenomeDir(object):
  """Create file indicating directory of each genome."""

  def __init__(self):
    pass

  def run(self, ncbi_genome_dir, user_genome_dir, db_prefix, output_file):
    """Create file indicating directory of each genome."""

    fout = open(output_file, 'w')

    for db in ['genbank', 'refseq']:
        for domain in ['archaea', 'bacteria']:
          domain_dir = os.path.join(ncbi_genome_dir, db, domain)
          for species_dir in os.listdir(domain_dir):
            full_species_dir = os.path.join(domain_dir, species_dir)
            if not os.path.isdir(full_species_dir):
                continue

            for assembly_dir in os.listdir(full_species_dir):
                accession = assembly_dir[0:assembly_dir.find('_', 4)]
                if db_prefix:
                    if db == 'genbank':
                        accession = 'GB_' + accession
                    else:
                        accession = 'RS_' + accession
            
                full_assembly_dir = os.path.join(full_species_dir, assembly_dir)
                fout.write('%s\t%s\n' % (accession, os.path.abspath(full_assembly_dir)))
                
    if user_genome_dir != 'None':
        for user in os.listdir(user_genome_dir):
            user_dir = os.path.join(user_genome_dir, user)
            if not os.path.isdir(user_dir):
                continue
                
            for genome_id in os.listdir(user_dir):
                full_path = os.path.join(user_dir, genome_id)
                fout.write('%s\t%s\n' % (genome_id, os.path.abspath(full_path)))

    fout.close()

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('ncbi_genome_dir', help='base directory leading to NCBI GenBank and RefSeq genomes.')
  parser.add_argument('user_genome_dir', help='base directory leading to user genomes; None to skip.')
  parser.add_argument('output_file', help='output metadata file')
  parser.add_argument('--db_prefix', help='append GTDB source prefixes to genome accessions', action="store_true")

  args = parser.parse_args()

  try:
    p = GenomeDir()
    p.run(args.ncbi_genome_dir, args.user_genome_dir, args.db_prefix, args.output_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
