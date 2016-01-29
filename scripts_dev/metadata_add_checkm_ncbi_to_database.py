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

__prog_name__ = 'metadata_add_checkm_ncbi_to_database.py'
__prog_desc__ = 'Add CheckM information for NCBI genomes to database.'

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

from biolib.taxonomy import Taxonomy

class AddCheckM(object):
  """Add CheckM information for NCBI genomes to database."""

  def __init__(self):
    self.metadata = {'Completeness': ['checkm_completeness', 'FLOAT'],
                      'Contamination': ['checkm_contamination', 'FLOAT'],
                      'Strain heterogeneity': ['checkm_strain_heterogeneity', 'FLOAT'],
                      'Marker lineage': ['checkm_marker_lineage', 'TEXT'],
                      '# genomes': ['checkm_genome_count', 'INT'],
                      '# markers': ['checkm_marker_count', 'INT'],
                      '# marker sets': ['checkm_marker_set_count', 'INT']}

  def run(self, checkm_qa_file):
    """Add CheckM data to database."""

    for header, data in self.metadata.iteritems():
      db_header, data_type = data

      temp_file = tempfile.NamedTemporaryFile(delete=False)
      with open(checkm_qa_file) as f:
        headers = f.readline().rstrip().split('\t')
        col_index = headers.index(header)

        for line in f:
          line_split = line.split('\t')
          genome_id = line_split[0]
          genome_id = genome_id[0:genome_id.find('_', 4)]

          if genome_id.startswith('GCA_'):
            genome_id = 'GB_' + genome_id
          elif genome_id.startswith('GCF_'):
            genome_id = 'RS_' + genome_id

          data = line_split[col_index]
          temp_file.write('%s\t%s\n' % (genome_id, data))

        temp_file.close()
        cmd = 'gtdb metadata import --table %s --field %s --type %s --metadatafile %s' % ('metadata_genes', db_header, data_type, temp_file.name)
        print cmd
        os.system(cmd)
        os.remove(temp_file.name)

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('checkm_qa_file', help='CheckM QA file for all genomes of interest')

  args = parser.parse_args()

  try:
    p = AddCheckM()
    p.run(args.checkm_qa_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
