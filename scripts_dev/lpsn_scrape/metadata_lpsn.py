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

__prog_name__ = 'metadata_lpsn.py'
__prog_desc__ = 'Produce metadata files describing type genera, species, and strains according to LPSN'

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
import ntpath
import argparse
from collections import defaultdict
import re

from biolib.taxonomy import Taxonomy


class MetadataLPSN(object):
    """Produce metadata files describing type genera, species, and strains according to LPSN."""

    def __init__(self):
        pass

    def run(self, lpsn_scrape_file, output_dir):
        """Create metadata by parsing assembly stats files."""

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # identify type genera, species, and strains according to LPSN
        fout_type_genera = open(os.path.join(output_dir, 'lpsn_genera.tsv'), 'w')
        fout_type_species = open(os.path.join(output_dir, 'lpsn_species.tsv'), 'w')
        fout_type_strains = open(os.path.join(output_dir, 'lpsn_strains.tsv'), 'w')

        fout_type_genera.write('lpsn_genus\tlpsn_type_genus\tlpsn_genus_authority\n')
        fout_type_species.write('lpsn_species\tlpsn_type_species\tlpsn_species_authority\n')
        fout_type_strains.write('lpsn_strain\n')

        strains = set()
        for line in open(lpsn_scrape_file):
            line_split = line.rstrip().split('\t')

            if line_split[0] == 'genus':
                genus = 'g__' + line_split[2]
                desc = line_split[3].strip()

                family = ''
                if int(line_split[1]) == 1:
                    family = desc[desc.find('family ?') + len('family ?'):].strip()
                    family = 'f__' + family.split()[0]

                desc = desc.replace(' ?', '')

                fout_type_genera.write('%s\t%s\t%s\n' % (genus, family, desc))
            elif line_split[0] == 'species':
                species = 's__' + line_split[2]
                desc = line_split[3].strip()

                genus = ''
                if int(line_split[1]) == 1:
                    genus = 'g__' + line_split[2].split()[0]

                fout_type_species.write('%s\t%s\t%s\n' % (species, genus, desc))
                processed_strains = []
                for i, strain in enumerate(line_split[4:]):
                    if i == 0 and strain.startswith('strain '):
                        strain = strain.replace("strain ", "", 1)
                    strain = re.sub(r'\(.+\)', ' ', strain)
                    strain = ' '.join(strain.split())
                    matchObj = re.match(r'^[\w|\s|\d|\.|-]+$', strain, re.M | re.I)
                    if matchObj:
                        processed_strains.append(strain)
                        if " " in strain:
                            strain = strain.replace(" ", "")
                            processed_strains.append(strain)
                fout_type_strains.write('{0} {1}\n'.format(line_split[2], "=".join(processed_strains)))


#         for strain_id in line_split[4:]:
#           strain = line_split[2] + ' ' + strain_id
#           if strain not in strains:
#             fout_type_strains.write('%s\n' % strain)
#             strains.add(strain)

        fout_type_genera.close()
        fout_type_species.close()
        fout_type_strains.close

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('lpsn_scrape_file', help='output file from scarping LPSN')
    parser.add_argument('output_dir', help='file containing LPSN metadata fot GTDB')

    args = parser.parse_args()

    try:
        p = MetadataLPSN()
        p.run(args.lpsn_scrape_file, args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
