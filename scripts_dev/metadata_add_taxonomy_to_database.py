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

__prog_name__ = 'metadata_add_taxonomy_to_database.py'
__prog_desc__ = 'Add genome tree taxonomy columns to the GTDB.'

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

class AddTaxonomy(object):
  """Add genome tre taxonomy to GTDB."""

  def __init__(self):
    pass

  def run(self, taxonomy_file, genome_list_file):
    """Add taxonomy to database."""

    genome_list = set()
    for line in open(genome_list_file):
      genome_list.add(line.rstrip().split('\t')[0])

    # read taxonomy file
    taxonomy = Taxonomy().read(taxonomy_file)

    # add full taxonomy string to database
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    for genome_id, taxa in taxonomy.iteritems():
      if genome_list_file and genome_id not in genome_list:
        continue

      taxa_str = ';'.join(taxa)
      temp_file.write('%s\t%s\n' % (genome_id, taxa_str))

    temp_file.close()
    cmd = 'gtdb metadata import --table %s --field %s --type %s --metadatafile %s' % ('metadata_taxonomy', 'gtdb_taxonomy', 'TEXT', temp_file.name)
    print cmd
    os.system(cmd)
    os.remove(temp_file.name)

    # add each taxonomic rank to database
    for i, rank in enumerate(Taxonomy.rank_labels):
      temp_file = tempfile.NamedTemporaryFile(delete=False)
      for genome_id, taxa in taxonomy.iteritems():
        if genome_list_file and genome_id not in genome_list:
          continue

        rank_str = taxa[i]
        if Taxonomy.rank_labels[i] == 'species':
          # ensure species name includes genus
          if taxa[i-1][3:] not in taxa[i]:
            rank_str = 's__' + taxa[i-1][3:] + ' ' + taxa[i][3:]

        temp_file.write('%s\t%s\n' % (genome_id, rank_str))

      temp_file.close()
      cmd = 'gtdb metadata import --table %s --field %s --type %s --metadatafile %s' % ('metadata_taxonomy', 'gtdb_' + rank, 'TEXT', temp_file.name)
      print cmd
      os.system(cmd)
      os.remove(temp_file.name)

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('taxonomy_file', help='taxonomy file with genome tree taxonomy strings')
  parser.add_argument('--genome_list', help='only process genomes in this list')

  args = parser.parse_args()

  try:
    p = AddTaxonomy()
    p.run(args.taxonomy_file, args.genome_list)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
