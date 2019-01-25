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

__prog_name__ = 'ncbi_genome_type.py'
__prog_desc__ = 'Identify genomes marked by NCBI as being a MAG or SAG.'

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
import argparse
from collections import defaultdict

class GenomeType(object):
    """Identify genomes marked by NCBI as being a MAG or SAG."""

    def __init__(self):
        pass

    def run(self, genome_file, output_file):
        """Identify genomes marked by NCBI as being a MAG or SAG."""
        
        fout = open(output_file, 'w')
        fout.write('genome_id\tncbi_genome_type\n')
        fout_sanity_check = open(output_file + '.raw', 'w')
        genome_count = 0
        for idx, line in enumerate(open(genome_file)):
            if idx % 100 == 0:
                sys.stdout.write('==> Processed %d genomes.\r' % idx)
                sys.stdout.flush()
                
            gid, genome_dir = line.strip().split('\t')
            
            if gid.startswith('U_'):
                fout.write('%s\t%s\n' % (gid, 'derived from metagenome'))
                continue

            assembly_id = os.path.basename(os.path.normpath(genome_dir))
            wgs_file = os.path.join(genome_dir, assembly_id + '_genomic.gbff')
            for line in open(wgs_file):
                if 'derived from metagenome' in line:
                    fout.write('%s\t%s\n' % (gid, 'derived from metagenome'))
                    fout_sanity_check.write(line)
                    break
                elif 'single cell' in line:
                    fout.write('%s\t%s\n' % (gid, 'derived from single cell'))
                    fout_sanity_check.write(line)
                    break
        sys.stdout.flush()
                    
        fout.close()
        fout_sanity_check.close()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_file', help='file indicating path to GTDB genome files')
    parser.add_argument('output_file', help='output file')

    args = parser.parse_args()

    try:
        p = GenomeType()
        p.run(args.genome_file, args.output_file)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
