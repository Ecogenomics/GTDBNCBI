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

__prog_name__ = 'ncbi_standardized_taxonomy.py'
__prog_desc__ = 'Parse NCBI taxonomy file to produce a standard 7-rank taxonomy for ammendable taxa.'

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
import traceback
from collections import namedtuple, defaultdict

from biolib.taxonomy import Taxonomy


class StandardizedTaxonomy(object):
  """Produce standardized 7-rank taxonomy file from NCBI taxonomy strings."""

  def __init__(self):
    pass

  def run(self, ncbi_taxonomy_file, output_consistent, output_inconsistent):
    """Produce standardized 7-rank taxonomy file from NCBI taxonomy strings."""

    fout_consistent = open(output_consistent, 'w')
    fout_inconsistent = open(output_inconsistent, 'w')
    for line in open(ncbi_taxonomy_file):
      line_split = line.strip().split('\t')

      assembly_accesssion = line_split[0]
      taxonomy = line_split[1].split(';')

      # remove unrecognized ranks (i.e., 'x__') and strain classification
      revised_taxonomy = []
      for t in taxonomy:
        if not t.startswith('x__') and not t.startswith('st__'):
          revised_taxonomy.append(t)

      # check if revised taxonomy is 7-ranks deep, has canonical
      # rank prefixes, and matched genus and species names
      if len(revised_taxonomy) == len(Taxonomy.rank_prefixes):
        valid_ranks = True
        for i, t in enumerate(revised_taxonomy):
          if not t.startswith(Taxonomy.rank_prefixes[i]):
            valid_ranks = False

        if valid_ranks:
          fout_consistent.write('%s\t%s\n' % (assembly_accesssion, ';'.join(revised_taxonomy)))
        else:
          fout_inconsistent.write('%s\t%s\n' % (assembly_accesssion, ';'.join(taxonomy)))
      else:
        fout_inconsistent.write('%s\t%s\n' % (assembly_accesssion, ';'.join(taxonomy)))

    fout_consistent.close()
    fout_inconsistent.close()

    print 'Genomes with a consistent taxonomy written to: %s' % output_consistent
    print 'Genomes with an inconsistent taxonomy written to: %s' % output_inconsistent

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('ncbi_taxonomy_file', help='NCBI tab-separated taxonomy file')
  parser.add_argument('output_consistent', help='output of genomes with a consistent taxonomy')
  parser.add_argument('output_inconsistent', help='output of genomes with an inconsistent taxonomy')

  args = parser.parse_args()

  try:
    p = StandardizedTaxonomy()
    p.run(args.ncbi_taxonomy_file, args.output_consistent, args.output_inconsistent)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
