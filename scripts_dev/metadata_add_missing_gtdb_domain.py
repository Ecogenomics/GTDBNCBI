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

__prog_name__ = 'metadata_add_missing_gtdb_domain.py'
__prog_desc__ = 'Fill empty gtdb_domain fields for NCBI genomes using the NCBI taxonomy.'

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
import csv
import argparse
import tempfile
from collections import defaultdict

from biolib.taxonomy import Taxonomy

csv.field_size_limit(sys.maxsize)

class FillMissingDomain(object):
  def __init__(self):
    pass

  def run(self, metadata_file):
    
    # find NCBI genomes with an empty gtdb_domain and valid NCBI taxonomy
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    
    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('genome')
            gtdb_domain_index = headers.index('gtdb_domain')
            ncbi_taxonomy_index = headers.index('ncbi_taxonomy')
            bHeader = False
        else:
            genome_id = row[genome_index]
            gtdb_domain = row[gtdb_domain_index]
            ncbi_taxonomy = row[ncbi_taxonomy_index]

            if (not gtdb_domain or gtdb_domain == 'd__') and ncbi_taxonomy:
                ncbi_domain = ncbi_taxonomy.split(';')[0]
                if ncbi_domain != 'd__':
                    print genome_id, ncbi_domain
                    temp_file.write('%s\t%s\n' % (genome_id, ncbi_domain))
                else:
                    print '[WARNING] NCBI genomes has no GTDB domain or valid NCBI taxonomy: %s' % genome_id

    temp_file.close()

    cmd = 'gtdb -r metadata import --table metadata_taxonomy --field gtdb_domain --type TEXT --metadatafile %s' % temp_file.name
    print cmd
    os.system(cmd)
    os.remove(temp_file.name)

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('metadata_file', help='metadata file for GTDB with information about all GTDB genomes')

  args = parser.parse_args()

  try:
    p = FillMissingDomain()
    p.run(args.metadata_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
