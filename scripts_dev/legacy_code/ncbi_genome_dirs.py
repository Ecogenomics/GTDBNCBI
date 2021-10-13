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
__version__ = '0.0.2'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse


class GenomeDir(object):
    """Create file indicating directory of each genome."""

    def __init__(self):
        pass

    def normaliseID(self, accession):
        normaccession = "G" + accession[4:accession.find('.', 0)]
        return normaccession

    def run(self, database_dir, output_file):
        """Create file indicating directory of each genome."""

        fout = open(output_file, 'w')

        for first_three in os.listdir(database_dir):
            onethird_species_dir = os.path.join(database_dir, first_three)
            # print onethird_species_dir
            if os.path.isfile(onethird_species_dir):
                continue
            for second_three in os.listdir(onethird_species_dir):
                twothird_species_dir = os.path.join(
                    onethird_species_dir, second_three)
                # print twothird_species_dir
                if os.path.isfile(twothird_species_dir):
                    continue
                for third_three in os.listdir(twothird_species_dir):
                    threethird_species_dir = os.path.join(
                        twothird_species_dir, third_three)
                    # print threethird_species_dir
                    if os.path.isfile(threethird_species_dir):
                        continue
                    for complete_name in os.listdir(threethird_species_dir):
                        full_path = os.path.join(
                            threethird_species_dir, complete_name)
                        if os.path.isfile(full_path):
                            continue
                        accession = complete_name[0:complete_name.find('_', 4)]

                        fout.write('{}\t{}\t{}\n'.format(accession, os.path.abspath(
                            full_path), self.normaliseID(accession)))

        fout.close()


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies')
    parser.add_argument('output_file', help='output metadata file')

    args = parser.parse_args()

    try:
        p = GenomeDir()
        p.run(args.genome_dir, args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
