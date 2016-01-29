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

__prog_name__ = 'taxonomy_fill.py'
__prog_desc__ = 'Ensure all taxonomy strings contain all 7 ranks.'

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

from biolib.taxonomy import Taxonomy 


class TaxonomyFill(object):

  def __init__(self):
    pass

  def run(self, input_taxonomy, output_taxonomy):
    fout = open(output_taxonomy, 'w')

    taxonomy = Taxonomy()
    t = taxonomy.read(input_taxonomy)
      
    for genome_id, taxon_list in t.iteritems():
        full_taxon_list = taxonmoy.fill_missing_ranks(taxon_list)
        fout.write('%s\t%s\n' % (genome_id, ';'.join(full_taxon_list)))
        
    fout.close()
          
if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input_taxonomy', help='input taxonomy file')
  parser.add_argument('output_taxonomy', help='output taxonomy file')

  args = parser.parse_args()

  try:
    p = TaxonomyFill()
    p.run(args.input_taxonomy, args.output_taxonomy)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
