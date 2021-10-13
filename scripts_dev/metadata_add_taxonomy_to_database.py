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
__version__ = '0.0.2'
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
    
  def truncate_taxonomy(self, metadata_file):
    """Truncate taxonomy string to just domain classification."""
    
    # get current GTDB taxonomy for all genomes
    gtdb_taxonomy = {}
    with open(metadata_file) as f:
        header = f.readline().strip().split('\t')
        
        gtdb_taxonomy_index = header.index('gtdb_taxonomy')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            gtdb_taxa = [t.strip() for t in line_split[gtdb_taxonomy_index].split(';')]
            gtdb_taxonomy[gid] = gtdb_taxa
            
    for i, rank in enumerate(Taxonomy.rank_labels):
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        for gid, taxa in list(gtdb_taxonomy.items()):
            if rank == 'domain':
                rank_str = taxa[i]
                temp_file.write('%s\t%s\n' % (gid, rank_str))
            else:
                temp_file.write('%s\t%s\n' % (gid, Taxonomy.rank_prefixes[i]))

        temp_file.close()
        cmd = 'gtdb -r metadata import --table %s --field %s --type %s --metadatafile %s' % ('metadata_taxonomy', 'gtdb_' + rank, 'TEXT', temp_file.name)
        print(cmd)
        os.system(cmd)
        os.remove(temp_file.name)
        
  def run(self, taxonomy_file, metadata_file, genome_list_file, truncate_taxonomy):
    """Add taxonomy to database."""
    
    if truncate_taxonomy:
        print('Truncating GTDB taxonomy to domain classification.')
        self.truncate_taxonomy(metadata_file)

    genome_list = set()
    if genome_list_file:
        for line in open(genome_list_file):
            if '\t' in line:
                genome_list.add(line.rstrip().split('\t')[0])
            else:
                genome_list.add(line.rstrip().split(',')[0])

    # read taxonomy file
    taxonomy = Taxonomy().read(taxonomy_file)

    # add each taxonomic rank to database
    for i, rank in enumerate(Taxonomy.rank_labels):
      temp_file = tempfile.NamedTemporaryFile(delete=False)
      for genome_id, taxa in list(taxonomy.items()):
        if genome_list_file and genome_id not in genome_list:
          continue

        rank_str = taxa[i]
        temp_file.write('%s\t%s\n' % (genome_id, rank_str))

      temp_file.close()
      cmd = 'gtdb -r metadata import --table %s --field %s --type %s --metadatafile %s' % ('metadata_taxonomy', 'gtdb_' + rank, 'TEXT', temp_file.name)
      print(cmd)
      os.system(cmd)
      os.remove(temp_file.name)

if __name__ == '__main__':
  print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
  print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('taxonomy_file', help='taxonomy file with genome tree taxonomy strings')
  parser.add_argument('metadata_file', help='metadata file for all genomes in GTDB')
  parser.add_argument('--genome_list', help='only process genomes in this list')
  parser.add_argument('--truncate_taxonomy', 
                        help='truncate current taxonomy strings to just the domain before updating taxonomy', 
                        action='store_true')

  args = parser.parse_args()

  try:
    p = AddTaxonomy()
    p.run(args.taxonomy_file,
            args.metadata_file,
            args.genome_list,
            args.truncate_taxonomy)
  except SystemExit:
    print("\nControlled exit resulting from an unrecoverable error or warning.")
  except:
    print("\nUnexpected error:", sys.exc_info()[0])
    raise
