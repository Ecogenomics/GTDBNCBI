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

__prog_name__ = 'ncbi_strain_checker.py'
__prog_desc__ = 'Check is all strain fields are consistent.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import os
import sys
import argparse
import tempfile
from collections import defaultdict
import re


class StrainChecker(object):
  """Extract genes in nucleotide space."""

  def __init__(self):
    pass

  def run(self,refsum,gensum,outdir):

    gbffreport_different = open(os.path.join(outdir,"gbffreport_different.txt"),'w')
    gbffreport_oneempty = open(os.path.join(outdir,"gbffreport_oneempty.txt"),'w')
    summaryreport = open(os.path.join(outdir,"strain_identifiers.txt"),'w')



    gbffreport_different.write("genome_id\tOrganism name\tncbi_strain_identifiers\tstrain (gbff)\tstrain (report)\tncbi_type_material\n")
    gbffreport_oneempty.write("genome_id\tOrganism name\tncbi_strain_identifiers\tstrain (gbff)\tstrain (report)\tncbi_type_material\n")
    summaryreport.write("genome_id\tOrganism name\tncbi_strain_identifiers\tstrain (gbff)\tstrain (report)\tncbi_type_material\n")
    
    list_file = [refsum,gensum]
    
    for sumfile in list_file:
        with open(sumfile,'r') as genomestrains :
            genomestrains.readline()
            for line in genomestrains:
                normal = True
                info = line.rstrip("\n").split("\t")
                filtered_list = filter(lambda a: a != '', [info[2],info[3],info[4]])
                if len(set(filtered_list)) > 1 :
                    gbffreport_different.write(line)
                    normal =False
                if info[2] == '' or info[2] == 'n/a':
                    info2 = "empty"
                else:
                    info2 = info[2]
                if info[3] == '' or info[3] == 'n/a':
                    info3 = "empty"
                else:
                    info3 = info[3]
                if info[4] == '' or info[4] == 'n/a':
                    info4 = "empty"
                else:
                    info4 = info[4]
                if len(set([info2,info3,info4])) > 1 and [info2,info3,info4].count('empty')>1 :
                    gbffreport_oneempty.write("{0}\n".format('\t'.join([info[0],info[1],info2,info3,info4])))
                    normal =False
                if normal:
                    summaryreport.write(line)
                
if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--refseq_summary', help='refseq summary generate from ncbi_strain_summary')
    parser.add_argument('--genbank_summary', help='genbank summary generate from ncbi_strain_summary')
    parser.add_argument('--out_dir', help='header that will be concatenated to the output files (gb,rs,gbk,....)')
  
    args = parser.parse_args()

    try:
        p = StrainChecker()
        p.run(args.refseq_summary,args.genbank_summary,args.out_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
            
            
        
        