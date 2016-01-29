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

__prog_name__ = 'ncbi_type_strain_metadata_table.py'
__prog_desc__ = 'Create table indicate which assemblies represent type strains.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import csv
import sys
import argparse
import traceback
from collections import namedtuple, defaultdict


class TypeStrainTable(object):
  """Create type strain metadata table."""

  def run(self, metadata_table, type_strain_file, output_file):
    """Create type strain metadata table."""

    # read type stain file
    type_strain_taxids = set()
    with open(type_strain_file) as f:
      f.readline()

      for line in f:
        line_split = line.split('\t')
        type_strain_taxid = line_split[0]
        type_strain_taxids.add(type_strain_taxid)

    # determine assemblies which represent type strains
    fout = open(output_file, 'w')
    fout.write('assembly_id\tncbi_type_strain\n')

    csv_reader = csv.reader(open(metadata_table))
    bHeader = True
    for row in csv_reader:
      if bHeader:
        taxid_index = row.index('ncbi_taxid')
        genome_index = row.index('genome')
        bHeader = False
      else:
        genome_id = row[genome_index]
        taxid = row[taxid_index]

        if taxid in type_strain_taxids:
          fout.write('%s\t%s\n' % (genome_id, 'yes'))

    fout.close()

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('gtdb_metadata_table', help='metadata table from GTDB indicating taxid for each assembly')
  parser.add_argument('type_strain_file', help='tab-separated values file indicating NCBI-derived type material')
  parser.add_argument('output_file', help='output file')

  args = parser.parse_args()

  try:
    p = TypeStrainTable()
    p.run(args.gtdb_metadata_table, args.type_strain_file, args.output_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
