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

__prog_name__ = 'create_ncbi_summary_file.py'
__prog_desc__ = 'Parse metadata file to generate a summary file. Similar to the assembly_summary.txt file in NCBI FTP'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import os
import sys
import argparse
import re

class AssemblyParser(object):
  """Identify assemblies exclusive to GenBank."""

  def __init__(self):
      pass

  def run(self, ncbi_directory_file, output_file , org_file):
      
      
    #===========================================================================
    # fout = open(output_file, 'w')
    # fout.write("# assembly_accession\ttaxid\n")
    # with open(ncbi_directory_file,'r') as dirfile:
    #     for line in dirfile:
    #         genomepath = line.split('\t')[1].strip()
    #         gid = line.split('\t')[0]
    #         items = os.listdir(genomepath)
    #         for names in items:
    #             if names.endswith("_assembly_report.txt"):
    #                 with open(os.path.join(genomepath,names)) as reportfile:
    #                     for reportline in reportfile:
    #                         if reportline.startswith("# Taxid:"):
    #                             reportline = reportline.replace("# Taxid:","").lstrip().rstrip()
    #                             fout.write("{0}\t{1}\n".format(gid,reportline))
    # fout.close()
    #===========================================================================
    
    foutorg = open(org_file, 'w')
    with open(ncbi_directory_file,'r') as dirfile:
        for line in dirfile:
            genomepath = line.split('\t')[1].strip()
            gid = line.split('\t')[0]
            prefix = 'GB_'
            if gid.startswith("GCF_"):
                prefix = 'RS_'
            items = os.listdir(genomepath)
            for names in items:
                if names.endswith("_assembly_report.txt"):
                    with open(os.path.join(genomepath,names)) as reportfile:
                        for reportline in reportfile:
                            if reportline.startswith("# Organism name:"):
                                reportline = reportline.replace("# Organism name:","").lstrip().rstrip()
                                reportline = re.sub(r'\([^)]*\)', '', reportline)
                                foutorg.write("{0}\t{1}\n".format(prefix+gid,reportline))
    foutorg.close()


if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('ncbi_directory_file', help='GenBank or Refseq directory file.')
  parser.add_argument('output_file', help='output file. (preferably in metadata folder)')
  parser.add_argument('org_file', help='output file listing the organism names. (preferably in metadata folder)')

  args = parser.parse_args()

  try:
    p = AssemblyParser()
    p.run(args.ncbi_directory_file, args.output_file, args.org_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
