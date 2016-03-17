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

__prog_name__ = 'ncbi_taxonomy.py'
__prog_desc__ = 'Parse NCBI taxonomy files to produce simplified summary file.'

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


class TaxonomyNCBI(object):
  """Parse NCBI taxonomy files to produce a simplified summary file."""

  def __init__(self):
    self.NodeRecord = namedtuple('NodeRecord', 'parent_tax_id rank division_id genetic_code_id')
    self.NameRecord = namedtuple('NamesRecord', 'name_txt')

    self.bacterial_division = '0'
    self.unassigned_division = '8'

  def _assembly_to_tax_id(self, assembly_metadata_file):
    """Determine taxonomic identifier for each assembly.

    Parameters
    ----------
    assembly_metadata_file : str
      Path to assembly metadata file.

    Returns
    -------
    dict : d[assembly_accession] -> tax_id
      Taxonomic identifier for each assembly.
    """

    d = {}
    with open(assembly_metadata_file) as f:
      headers = f.readline().split('\t')

      taxid_index = headers.index('taxid')

      for line in f:
        line_split = line.split('\t')
        assembly_accession = line_split[0]
        taxid = line_split[taxid_index]

        if assembly_accession in d:
          print '[Error] Duplicate assembly accession: %s' % assembly_accession
          sys.exit(-1)

        d[assembly_accession] = taxid

    return d

  def _read_nodes(self, nodes_file):
    """Read NCBI nodes.dmp file.

    Parameters
    ----------
    nodes_file : str
      Path to NCBI nodes.dmp file.

    Returns
    -------
    dict : d[tax_id] -> NodeRecord
      Node record for all nodes.
    """

    d = {}
    for line in open(nodes_file):
      line_split = [t.strip() for t in line.split('|')]

      tax_id = line_split[0]
      parent_tax_id = line_split[1]
      rank = line_split[2]
      division_id = line_split[4]
      genetic_code_id = line_split[6]

      d[tax_id] = self.NodeRecord(parent_tax_id, rank, division_id, genetic_code_id)

    return d

  def _read_names(self, names_file):
    """Read NCBI names.dmp file.

    Parameters
    ----------
    names_file : str
      Path to NCBI names.dmp file.

    Returns
    -------
    dict : d[tax_id] -> NameRecord
      Name record of nodes marked as 'scientific name'.
    """

    d = {}
    for line in open(names_file):
      line_split = [t.strip() for t in line.split('|')]

      tax_id = line_split[0]
      name_txt = line_split[1]
      unique_name = line_split[2]
      name_class = line_split[3]

      if name_class == 'scientific name':
        d[tax_id] = self.NameRecord(name_txt)

    return d

  def run(self, taxonomy_dir, assembly_metadata_file, output_file):
    """Read NCBI taxonomy information and create summary output files."""

    # parse metadata file and taxonomy files
    assembly_to_tax_id = self._assembly_to_tax_id(assembly_metadata_file)

    node_records = self._read_nodes(os.path.join(taxonomy_dir, 'nodes.dmp'))
    print 'Read %d node records.' % len(node_records)

    name_records = self._read_names(os.path.join(taxonomy_dir, 'names.dmp'))
    print 'Read %d name records.' % len(name_records)

    # traverse taxonomy tree for each assembly
    fout = open(output_file, 'w')

    print 'Number of assemblies: %d' % len(assembly_to_tax_id)
    for assembly_accession, tax_id in assembly_to_tax_id.iteritems():
      # traverse taxonomy tree to the root which is 'cellular organism' for genomes,
      # 'other suquences' for plasmids, and 'unclassified sequences' for metagenomic libraries
      taxonomy = []
      cur_tax_id = tax_id

      if cur_tax_id not in name_records:
        print '[Warning] Assembly %s has an invalid taxid: %s' % (assembly_accession, tax_id)
        continue

      roots = ['cellular organisms', 'other sequences', 'unclassified sequences', 'Viruses', 'Viroids']
      while name_records[cur_tax_id].name_txt not in roots:
        if cur_tax_id == '1':
          print '[Error] TaxId %s reached root of taxonomy tree: %s' % (tax_id, taxonomy)
          sys.exit(-1)

        try:
          node_record = node_records[cur_tax_id]

          if node_record.rank in Taxonomy.rank_labels:
            rank_index = Taxonomy.rank_labels.index(node_record.rank)
            rank_prefix = Taxonomy.rank_prefixes[rank_index]
          else:
            # unrecognized rank
            rank_prefix = 'x__'
            if node_record.rank == 'superkingdom':
              rank_prefix = 'd__'
            elif node_record.rank == 'no rank' and len(taxonomy) == 0:
              # the taxonomy id is likely for a specific strain,
              # so mark it as a strain
              rank_prefix = 'st__'

          taxonomy.append(rank_prefix + name_records[cur_tax_id].name_txt)

          cur_tax_id = node_record.parent_tax_id
        except:
          print traceback.format_exc()
          print taxonomy

      taxonomy.reverse()
      taxa_str = ';'.join(taxonomy)
      fout.write('%s\t%s\n' % (assembly_accession, taxa_str))

    fout.close()

    print ''
    print 'Taxonomy information for assemblies written to: %s' % output_file

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('taxonomy_dir', help='directory containing NCBI taxonomy files')
  parser.add_argument('assembly_metadata_file', help='file with metadata for each assembly obtain from NCBI FTP site')
  parser.add_argument('output_file', help='output taxonomy file')

  args = parser.parse_args()

  try:
    p = TaxonomyNCBI()
    p.run(args.taxonomy_dir, args.assembly_metadata_file, args.output_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
