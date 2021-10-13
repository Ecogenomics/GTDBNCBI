#!/usr/bin/python


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

import os
import sys
import argparse

from biolib.seq_io import read_fasta

__prog_name__ = 'generate_constraint_file_for_fasttree.py'
__prog_desc__ = 'Generate a constraint File for Fasttree .'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au.au'
__status__ = 'Development'


class Constraintgenerator(object):
    def __init__(self):
        pass

    def run(self,msa_file,constraint_dir,outfile):
        msa_dict = read_fasta(msa_file)
        outdict = dict((key, []) for key in msa_dict.keys())
        onlyfiles = [os.path.join(constraint_dir, f) for f in os.listdir(constraint_dir) if os.path.isfile(os.path.join(constraint_dir, f))]
        for constraintfile in onlyfiles:
            constraintlist = []
            with open(constraintfile) as f:
                for line in f:
                    constraintlist.append(line.strip())
                for k,v in outdict.items():
                    if k in constraintlist:
                        outdict.get(k).append('1')
                    else:
                        outdict.get(k).append('0')
        outf = open(outfile,'w')
        for outk,outval in outdict.items():
            outf.write(">{0}\n{1}\n".format(outk,''.join(outval)))
        outf.close()


if __name__ == '__main__':
  print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
  print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--msa_file', help = 'msa file ' )
  parser.add_argument('--constraint_dir', help='directory of all files listing a constraint for Fasttree.')
  parser.add_argument('--outfile', help='Output file.')


  args = parser.parse_args()

  try:
    p = Constraintgenerator()
    p.run(args.msa_file,args.constraint_dir,args.outfile)
  except SystemExit:
    print("\nControlled exit resulting from an unrecoverable error or warning.")
  except:
    print("\nUnexpected error:", sys.exc_info()[0])
    raise
