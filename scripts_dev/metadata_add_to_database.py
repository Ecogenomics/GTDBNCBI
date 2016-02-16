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

__prog_name__ = 'metadata_add_to_database.py'
__prog_desc__ = 'Add columns in a metadata file to the GTDB.'

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
import tempfile
from collections import defaultdict

class AddMetadata(object):
  """Add all columns in a metadata file to the GTDB."""

  def __init__(self):
    pass

  def run(self, metadata_file, metadata_desc_file, genome_list_file):
    """Add metadata."""

    # get genomes to process
    genome_list = set()
    if genome_list_file:
      for line in genome_list_file:
        genome_list.add(line.rstrip().split('\t')[0])

    # get database table and data type of each metadata field
    metadata_type = {}
    metadata_table = {}
    for line in open(metadata_desc_file):
      line_split = line.strip().split('\t')
      metadata_type[line_split[0]] = line_split[2]
      metadata_table[line_split[0]] = line_split[3]

    # read metadata file
    metadata = defaultdict(lambda : defaultdict(str))
    with open(metadata_file) as f:
      fields = [x.strip() for x in f.readline().split('\t')]

      for line in f:
        line_split = line.rstrip().split('\t')

        genome_id = line_split[0]
        for i, value in enumerate(line_split[1:]):
          metadata[fields[i+1]][genome_id] = value

    print metadata.keys()

    # add each field to the database
    for field in metadata:
      temp_file = tempfile.NamedTemporaryFile(delete=False)

      if field not in metadata_type:
        continue

      data_type = metadata_type[field]
      table = metadata_table[field]

      for genome_id, value in metadata[field].iteritems():
        if genome_id == 'GCF_000826165.1':
          continue

        try:
          if float(value) and data_type in ['INT', 'INTEGER']:
            # assume specified data type is correct and that we may need
            # to cast floats to integers
            value = str(int(float(value)))
        except:
          pass

        if value.strip():
          if genome_id.startswith('GCA_'):
            genome_id = 'GB_' + genome_id
          elif genome_id.startswith('GCF_'):
            genome_id = 'RS_' + genome_id

          if not genome_list or genome_id in genome_list:
            temp_file.write('%s\t%s\n' % (genome_id, value))

      temp_file.close()
      cmd = 'gtdb -r metadata import --table %s --field %s --type %s --metadatafile %s' % (table, field, data_type, temp_file.name)
      print cmd
      os.system(cmd)
      os.remove(temp_file.name)

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('metadata_file', help='tab-separated values table with metadata')
  parser.add_argument('metadata_desc_file', help='tab-separated values table with description of metadata fields')
  parser.add_argument('--genome_list', help='only process genomes in this list', default=None)

  args = parser.parse_args()

  try:
    p = AddMetadata()
    p.run(args.metadata_file, args.metadata_desc_file, args.genome_list)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
