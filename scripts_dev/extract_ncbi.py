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

__prog_name__ = 'extract_ncbi.py'
__prog_desc__ = 'Extract NCBI genome information.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse


class Extract(object):
  """Extract NCBI genome information.

  This script assumes that genomes are stored in individual
  directories in the following format:

  <domain>/<genome_dir>/latest_assembly_versions/<assembly_id>/<assembly_id>.<file_type>.gz

  where <domain> is archaea or bacteria;

  This is the directory structure used by NCBI.
  """

  def __init__(self):
    pass

  def run(self, ncbi_dir):
    # check if NCBI directory looks correct
    dirs = os.listdir(ncbi_dir)
    if 'bacteria' not in dirs or 'archaea' not in dirs:
      print "[Error] Expect 'bacteria' and 'archaea' directories."
      sys.exit()

    for domain in ['archaea']: #, 'bacteria']:

      domain_dir = os.path.join(ncbi_dir, domain)
      for d in os.listdir(domain_dir):

        genome_dir = os.path.join(domain_dir, d)
        if os.path.isdir(genome_dir):
          assembly_dirs = os.listdir(genome_dir)
          for assembly_id in assembly_dirs:
            assembly_dir = os.path.join(genome_dir, assembly_id)
            os.system('gunzip %s' % os.path.join(assembly_dir, '*.gz'))

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('ncbi_dir', help='directory with NCBI genomes')

  args = parser.parse_args()

  try:
    p = Extract()
    p.run(args.ncbi_dir)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
