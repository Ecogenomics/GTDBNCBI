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

__prog_name__ = 'ncbi_unique_assemblies.py'
__prog_desc__ = 'Count number of unique Archaea and Bacterial assemblies.'

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

class UniqueAssemblies(object):
  """Count number of unique Archaea and Bacterial assemblies."""

  def __init__(self):
    pass

  def run(self, genome_dir, output_file):
    """Create metadata by parsing assembly stats files."""

    fout = open(output_file, 'w')

    globally_unique = set()
    for domain in ['archaea', 'bacteria']:
      assemblies = defaultdict(list)

      domain_dir = os.path.join(genome_dir, domain)
      for species_dir in os.listdir(domain_dir):
        full_species_dir = os.path.join(domain_dir, species_dir)
        for assembly_dir in os.listdir(full_species_dir):
          accession = assembly_dir[0:assembly_dir.find('_', 4)]
          assemblies[accession].append(species_dir)

          if accession in globally_unique:
            print 'problem', accession
          globally_unique.add(accession)

      print 'Identified %d unique %s assemblies' % (len(assemblies), domain.capitalize())

      for accession, species_dirs in assemblies.iteritems():
        fout.write('%s\t%s\t%s\n' % (domain.capitalize(), accession, ','.join(species_dirs)))

    fout.close()

    print 'Identified %d globally unique assembly accession numbers.' % len(globally_unique)

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies')
  parser.add_argument('output_file', help='output file of all unique assemblies')

  args = parser.parse_args()

  try:
    p = UniqueAssemblies()
    p.run(args.genome_dir, args.output_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
