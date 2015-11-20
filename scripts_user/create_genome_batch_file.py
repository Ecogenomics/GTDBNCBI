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

__prog_name__ = 'create_genome_batch_file'
__prog_desc__ = 'Create genome batch file for GTDB.'

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

class GenomeBatchFile(object):
    def __init__(self):
        pass

    def run(self, binFolder, studyDesc, outFile, extension):
      fout = open(outFile, 'w')

      files = os.listdir(binFolder)

      binCount = 0
      for f in files:
        if f.endswith(extension):
          binId = f[0:f.rfind('.')]
          fout.write('%s\t%s\t%s\n' % (os.path.join(binFolder, f), binId, binId + ' (' + studyDesc + ')'))
          binCount += 1

      if not binCount:
        print 'No bins identified in %s. Check that your files have the correct extension (-x).' % binFolder
      else:
        print 'Process %d bins.' % binCount

      fout.close()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bin_folder', help='directory containing bins (please give your bins informative names when possible)')
    parser.add_argument('study_desc', help='brief description of your study/sample/project to associate with all bins')
    parser.add_argument('out_file', help='output file')
    parser.add_argument('-x', '--extension', help='extension of bins', default = 'fna')

    args = parser.parse_args()

    try:
        genomeBatchFile = GenomeBatchFile()
        genomeBatchFile.run(args.bin_folder, args.study_desc, args.out_file, args.extension)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
