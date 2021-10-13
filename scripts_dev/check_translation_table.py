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

__prog_name__ = 'check_translation_table.py'
__prog_desc__ = 'Check translation table used by Prodigal.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2018'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import shutil
import ntpath
import argparse


class Check(object):
    """Check translation table used by Prodigal.""" 

    def __init__(self):
        pass

    def run(self, gtdb_genome_dir, output_file):
        """Check translation table used by Prodigal."""
        
        fout = open(output_file, 'w')

        # get path to all unprocessed genome files
        for domain in ['archaea', 'bacteria']:
            print('Checking domain: %s' % domain)
            domain_dir = os.path.join(gtdb_genome_dir, domain)
            
            for idx, sp in enumerate(os.listdir(domain_dir)):
                print(idx, sp)
                sp_dir = os.path.join(domain_dir, sp)
                
                for genome in os.listdir(sp_dir):
                    genome_dir = os.path.join(sp_dir, genome)

                    try:
                        gff = os.path.join(genome_dir, genome + '_genomic.gff')
                        if os.path.exists(gff):
                            ncbi_transl_table = None
                            for line in open(gff):
                                if line[0] == '#':
                                    continue
                                    
                                if 'transl_table=' in line:
                                    trans_num = line[line.rfind('=')+1:].strip()
                                    ncbi_transl_table = int(trans_num)
                                    break
                            
                            if ncbi_transl_table: # not all GenBank genomes have the translation table indicated
                                prodigal_file = os.path.join(genome_dir, 'prodigal', 'prodigal_translation_table.tsv')
                                prodigal_transl_table = int(open(prodigal_file).readline().strip().split('\t')[1])
                                
                                if ncbi_transl_table != prodigal_transl_table:
                                    fout.write('%s\t%d\t%d\n' % (genome, ncbi_transl_table, prodigal_transl_table))
                                    print('%s\t%d\t%d' % (genome, ncbi_transl_table, prodigal_transl_table))
                    except:
                        print("\nUnexpected error:", sys.exc_info()[0])
                        print('Failed on: %s' % genome_dir)
                        sys.exit(-1)

        fout.close()

if __name__ == '__main__':
  print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
  print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('gtdb_genome_dir', help='directory leading to archaeal and bacterial genomes')
  parser.add_argument('output_file', help='output file')

  args = parser.parse_args()

  try:
    p = Check()
    p.run(args.gtdb_genome_dir, args.output_file)
  except SystemExit:
    print("\nControlled exit resulting from an unrecoverable error or warning.")
  except:
    print("\nUnexpected error:", sys.exc_info()[0])
    raise
