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

__prog_name__ = 'ncbi_genome_dirs.py'
__prog_desc__ = 'Produce file indicating the directory of each genome.'

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
from collections import defaultdict


class GenomeDir(object):
    """Create file indicating directory of each genome."""

    def __init__(self):
        pass

    def run(self, genome_dir, output_file, ftp_dir):
        """Create file indicating directory of each genome."""

        fout = open(output_file, 'w')

        for domain in ['archaea', 'bacteria']:
            domain_dir = os.path.join(genome_dir, domain)
            for species_dir in os.listdir(domain_dir):
                full_species_dir = os.path.join(domain_dir, species_dir)
                if os.path.isfile(full_species_dir):
                    continue
                if not os.listdir(full_species_dir):
                    continue
                if ftp_dir:
                    full_species_dir = os.path.join(full_species_dir, "latest_assembly_versions")
                if os.path.isdir(full_species_dir):
                    for assembly_dir in os.listdir(full_species_dir):
                        accession = assembly_dir[0:assembly_dir.find('_', 4)]
                        full_assembly_dir = os.path.join(full_species_dir, assembly_dir)
                        fout.write('%s\t%s\n' % (accession, os.path.abspath(full_assembly_dir)))

        fout.close()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies')
    parser.add_argument('output_file', help='output metadata file')
    parser.add_argument('-ftp', dest='ftp_directories', action='store_true',
                        help='Use if the FTP directories are parsed ( latest_assembly_versions')

    args = parser.parse_args()

    try:
        p = GenomeDir()
        p.run(args.genome_dir, args.output_file, args.ftp_directories)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
