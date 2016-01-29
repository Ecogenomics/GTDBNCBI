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

__prog_name__ = 'ncbi_type_strains.py'
__prog_desc__ = 'Determine type strains.'

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

class NameRecord(object):
  def __init__(self):
    self.name_txt = None
    self.type_material = set()

class TypeStrains(object):
  """Determine type strains.

    Type strains are determined by identifying the 'type material'
    of species identifiers in the NCBI taxonomy.
  """

  def __init__(self):
    self.NodeRecord = namedtuple('NodeRecord', 'parent_tax_id rank division_id genetic_code_id')

    self.bacterial_division = '0'
    self.unassigned_division = '8'

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

    d = defaultdict(lambda : NameRecord())
    for line in open(names_file):
      line_split = [t.strip() for t in line.split('|')]

      tax_id = line_split[0]
      name_txt = line_split[1]
      unique_name = line_split[2]
      name_class = line_split[3]

      if name_class == 'scientific name':
        d[tax_id].name_txt = name_txt
      elif name_class == 'type material':
        d[tax_id].type_material.add(name_txt)

    return d

  def run(self, taxonomy_dir, output_file):
    """Determine type strains."""

    node_records = self._read_nodes(os.path.join(taxonomy_dir, 'nodes.dmp'))
    print 'Read %d node records.' % len(node_records)

    name_records = self._read_names(os.path.join(taxonomy_dir, 'names.dmp'))
    print 'Read %d name records.' % len(name_records)

    # identify type strains
    type_strain_tax_ids = set()
    for tax_id, node_record in node_records.iteritems():
      parent_node_record = node_records[node_record.parent_tax_id]

      if parent_node_record.rank == 'species':
        parent_name_record = name_records[node_record.parent_tax_id]
        name_record = name_records[tax_id]

        for type_material in parent_name_record.type_material:
          if type_material in name_record.name_txt:
            type_strain_tax_ids.add(tax_id)

    fout = open(output_file, 'w')
    fout.write('Strain taxid\tStrain name\tSpecies taxid\tSpecies name\n')
    for tax_id in type_strain_tax_ids:
      node_record = node_records[tax_id]
      name_record = name_records[tax_id]
      parent_name_record = name_records[node_record.parent_tax_id]
      fout.write('%s\t%s\t%s\t%s\n' % (tax_id, name_record.name_txt, node_record.parent_tax_id, parent_name_record.name_txt))
    fout.close()

    print ''
    print 'The %d identified type strains written to: %s' % (len(type_strain_tax_ids), output_file)

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('taxonomy_dir', help='directory containing NCBI taxonomy files')
  parser.add_argument('output_file', help='output taxonomy file')

  args = parser.parse_args()

  try:
    p = TypeStrains()
    p.run(args.taxonomy_dir, args.output_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
