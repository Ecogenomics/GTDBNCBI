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

__prog_name__ = 'metadata_add_predicted_domain_to_database.py'
__prog_desc__ = 'Add predicted domain of genomes to the GTDB.'

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


class AddDomain(object):
  """Add predicted domain to GTDB."""

  def __init__(self):
    pass

  def run(self, domain_report):
    """Add predicted domain to GTDB."""
    
    domains = {}
    with open(domain_report) as f:
        f.readline()
        
        for line in f:
            line_split = line.strip().split('\t')
            genome_id = line_split[0]
            predicted_domain = line_split[1]
            gtdb_taxonomy = line_split[5]
            gtdb_domain = gtdb_taxonomy.split(';')[0]

            if predicted_domain != 'None':
                if gtdb_domain == 'd__' or predicted_domain == gtdb_domain:
                    domains[genome_id] = predicted_domain

    # add predicted domains to gtdb_domain field
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    for genome_id, domain in list(domains.items()):
        temp_file.write('%s\t%s\n' % (genome_id, domain))
    temp_file.close()
    
    cmd = 'gtdb -r metadata import --table metadata_taxonomy --field gtdb_domain --type TEXT --metadatafile %s' % temp_file.name
    print(cmd)
    os.system(cmd)
    #os.remove(temp_file.name)

if __name__ == '__main__':
  print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
  print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('domain_report', help='taxonomy file with genome tree taxonomy strings')

  args = parser.parse_args()

  try:
    p = AddDomain()
    p.run(args.domain_report)
  except SystemExit:
    print("\nControlled exit resulting from an unrecoverable error or warning.")
  except:
    print("\nUnexpected error:", sys.exc_info()[0])
    raise
