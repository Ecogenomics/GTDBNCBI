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
__prog_desc__ = 'Parse NCBI taxonomy file to produce a standard 7-rank taxonomy for amenable taxa.'

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

        # create longest taxonomy string possible with canonical ranks
        canonical_taxonomy = {}
        for i, taxon in enumerate(revised_taxonomy):
            rank_prefix = taxon[0:3]
            if rank_prefix in Taxonomy.rank_prefixes:
                if rank_prefix == 's__':
                    valid_name, reason = Taxonomy().validate_species_name(taxon)
                    if valid_name:
                        canonical_taxonomy[Taxonomy.rank_prefixes.index(rank_prefix)] = taxon
                else:
                    canonical_taxonomy[Taxonomy.rank_prefixes.index(rank_prefix)] = taxon
                
        # fill in missing ranks where possible
        if canonical_taxonomy:
            for i in xrange(0, max(canonical_taxonomy.keys())):
                if i in canonical_taxonomy and (i+1) not in canonical_taxonomy:
                    taxon = canonical_taxonomy[i][3:]
                    if taxon[0] == '{':
                        canonical_taxonomy[i+1] = Taxonomy.rank_prefixes[i+1] + taxon.replace(Taxonomy.rank_labels[i]+'}', 
                                                                                                Taxonomy.rank_labels[i+1]+'}')
                    else:
                        canonical_taxonomy[i+1] = Taxonomy.rank_prefixes[i+1] + '{undefined %s %s}' % (taxon, Taxonomy.rank_labels[i+1])

        cur_taxonomy = []
        for i in xrange(0, len(Taxonomy.rank_prefixes)):
            if i in canonical_taxonomy:
                cur_taxonomy.append(canonical_taxonomy[i])
            else:
                break # unable to correctly determine a valid taxonomy below this rank

        if len(cur_taxonomy) > 0:
            if len(cur_taxonomy) != len(Taxonomy.rank_prefixes):
                cur_taxonomy = cur_taxonomy + list(Taxonomy.rank_prefixes[len(cur_taxonomy):])
            fout_consistent.write('%s\t%s\n' % (assembly_accesssion, ';'.join(cur_taxonomy)))
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
