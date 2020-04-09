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

__prog_name__ = 'genomes_to_process_metadata.py'
__prog_desc__ = 'Determine new or modified genomes requiring processing.'

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


class GenomeToProcess(object):
  """Determine new or modified genomes requiring processing.

  The status of each genome is determined by comparing the
  latest rsync from NCBI with the existing genomes. A detailed
  change long is created indicating all modification. This
  script determines which genomes require metadata to be
  recalculated as they are:
    1. a new genome
    2. a modified genomes with new nucleotide data
    3. a genome with new gene calling
  """

  def __init__(self):
    pass

  def run(self, change_log, output_file):
    fout = open(output_file, 'w')

    assembles_to_process = {}
    for line in open(change_log):
      line_split = line.rstrip().split(':')

      if line_split[0] == 'CREATE':



if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('change_log', help='change log produced by comparing rsync directory with current genomes')
  parser.add_argument('output_file', help='output file with path to all genomes requiring metadata to be calculated')

  args = parser.parse_args()

  try:
    p = GenomeToProcess()
    p.run(args.change_log, args.output_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
