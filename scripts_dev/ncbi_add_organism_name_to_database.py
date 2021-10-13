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

__prog_name__ = 'ncbi_add_organism_name_to_database.py'
__prog_desc__ = 'Add NCBI organism name to the GTDB.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2017'
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


class AddOrganismName(object):
  """Add organism name to GTDB."""

  def __init__(self):
    pass

  def run(self, organism_name_file, genome_list_file):
    """Add organism name to database."""

    genome_list = set()
    if genome_list_file:
        for line in open(genome_list_file):
            if '\t' in line:
                genome_list.add(line.rstrip().split('\t')[0])
            else:
                genome_list.add(line.rstrip().split(',')[0])

    # add full taxonomy string to database
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    for line in open(organism_name_file):
        line_split = line.strip().split('\t')
        
        gid = line_split[0]
        org_name = line_split[1]
        if genome_list_file and gid not in genome_list:
            continue

        temp_file.write('%s\t%s\n' % (gid, org_name))

    temp_file.close()
    cmd = 'gtdb -r metadata import --table %s --field %s --type %s --metadatafile %s' % ('metadata_ncbi', 'ncbi_organism_name', 'TEXT', temp_file.name)
    print(cmd)
    os.system(cmd)
    os.remove(temp_file.name)

if __name__ == '__main__':
  print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
  print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('organism_name_file', help='file with organism name for each genome')
  parser.add_argument('--genome_list', help='only process genomes in this list')

  args = parser.parse_args()

  try:
    p = AddOrganismName()
    p.run(args.organism_name_file, args.genome_list)
  except SystemExit:
    print("\nControlled exit resulting from an unrecoverable error or warning.")
  except:
    print("\nUnexpected error:", sys.exc_info()[0])
    raise
