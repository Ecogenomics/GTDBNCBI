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

import os
import sys
import re
import argparse
from collections import defaultdict
import random
import string

__prog_name__ = 'extract_metadata_from_gbff.py'
__prog_desc__ = 'Produce metadata file from NCBI gbff file.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2016'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@qfab.org'
__status__ = 'Development'


class Metadata(object):
    """Create metadata file from the genome.gbff file for each NCBI assembly."""

    def __init__(self):
        self.gbff_fields = ['isolation_source', 'country', 'lat_lon']

        self.stats_info = {}

    def _randomword(self, length):
        return ''.join(random.choice(string.lowercase) for i in range(length))

    def _parse_gbff(self, gbff_file):
        metadata_gbff = [''] * len(self.gbff_fields)

        if not os.path.exists(gbff_file):
            return metadata_gbff

        pattern_gene = re.compile("^\s{0,20}\w")
        pattern_source = re.compile("^\s{5}source\s{10}")
        source_info_bool = False
        randomstring = self._randomword(10)
        source_info = []
        # We read the file line by line
        for line in open(gbff_file, 'r'):
            # when we reach the source part of the genbank file,we start
            # storing the line
            if pattern_source.match(line):
                source_info_bool = True
            # when the source paragraph is finished we stop reading the file
            elif pattern_gene.match(line) and source_info_bool:
                break
            # we store each linein an array
            elif source_info_bool:
                # we replace all the / character by a random string except the first one
                # this is because the / is used to separate the different metadata ex /organism, /lat_lon
                #/ will be used to separate those metadata later on

                line = re.sub(
                    r"(?!^\/)\/", randomstring, ' '.join(line.split()))
                source_info.append("{0} ".format(' '.join(line.split())))
        # all the lines from source are then concatenated in one and re split
        # following the /
        source_info_string = ''.join(x for x in source_info)
        # we redeclare source_info array
        source_info = source_info_string.split("/")
        source_info_dict = {}
        # all the metadata should be separated in each cell of the array.
        # we construct a dictionary from the array
        for info in source_info:
            if "=" in info:
                try:
                    k, v = info.split("=", 1)
                except:
                    print info
                source_info_dict[k] = v

        for field in self.gbff_fields:
            if field in source_info_dict:
                # we store the field replacing the special characters that can disturb the databases or the csv output
                # we replace the random string by /
                metadata_gbff[self.gbff_fields.index(field)] = source_info_dict[field].replace(
                    '"', '').replace(',', ';').replace("'", " ").replace(randomstring, "/").rstrip()
        return metadata_gbff

    def run(self, genome_dir, output_file):
        """Create metadata by parsing assembly gbff files."""

        fout = open(output_file, 'w')
        fout.write('Assembly accession')
        fout.write(
            '\t' + '\t'.join(['ncbi_' + x.lower() for x in self.gbff_fields]))
        fout.write('\n')

        processed_assemblies = defaultdict(list)
        count = 1
        for domain in ['archaea', 'bacteria']:
            domain_dir = os.path.join(genome_dir, domain)
            for species_dir in os.listdir(domain_dir):
                full_species_dir = os.path.join(domain_dir, species_dir)
                for assembly_dir in os.listdir(full_species_dir):
                    print count
                    count += 1
                    accession = assembly_dir[0:assembly_dir.find('_', 4)]

                    processed_assemblies[accession].append(species_dir)
                    if len(processed_assemblies[accession]) >= 2:
                        print '%s\t%s' % (accession, ','.join(processed_assemblies[accession]))
                        continue

                    full_assembly_dir = os.path.join(
                        full_species_dir, assembly_dir)

                    assembly_stat_file = os.path.join(
                        full_assembly_dir, assembly_dir + '_assembly_stats.txt')
                    if os.path.exists(assembly_stat_file):

                        gbff_file = os.path.join(
                            full_assembly_dir, assembly_dir + '_genomic.gbff')
                        gbff_stats = self._parse_gbff(gbff_file)
                        fout.write(accession + '\t%s' %
                                   ('\t'.join(map(str, gbff_stats))))
                    fout.write('\n')

        fout.close()


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies')
    parser.add_argument('output_file', help='output metadata file')

    args = parser.parse_args()

    try:
        p = Metadata()
        p.run(args.genome_dir, args.output_file)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
