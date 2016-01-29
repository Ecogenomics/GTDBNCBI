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

__prog_name__ = 'genbank_only_assemblies.py'
__prog_desc__ = 'Identify assemblies in GenBank that are not in RefSeq.'

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


class GenBankExclusive(object):
  """Identify assemblies exclusive to GenBank."""

  def __init__(self):
    pass

  def run(self, gb_bac_assembly_file, gb_ar_assembly_file, output_file):

    # identify assemblies exclusive to GenBank
    print 'Identifying assemblies exclusive to GenBank.'
    fout = open(output_file, 'w')
    count = 0
    for assembly_file in (gb_bac_assembly_file, gb_ar_assembly_file):
      with open(assembly_file) as f:
        headers = f.readline().rstrip().split('\t')
        gbrs_paired_asm_index = headers.index('gbrs_paired_asm')

        for line in f:
          line_split = line.rstrip().split('\t')
          gbrs_pair = line_split[gbrs_paired_asm_index]
          if gbrs_pair.strip() == 'na':
            fout.write(line_split[0] + '\n')
            count += 1

    print '  Identified %d assemblies.' % count

    fout.close()

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('gb_bac_assembly_file', help='GenBank bacterial assemblies metadata file from NCBI FTP site')
  parser.add_argument('gb_ar_assembly_file', help='GenBank archaeal assemblies metadata file from NCBI FTP site')
  parser.add_argument('output_file', help='output file with GenBank exclusive assemblies')

  args = parser.parse_args()

  try:
    p = GenBankExclusive()
    p.run(args.gb_bac_assembly_file, args.gb_ar_assembly_file, args.output_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
