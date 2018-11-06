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

__prog_name__ = 'ncbi_strain_summary.py'
__prog_desc__ = 'Parse the assembly report file, the genomic.gbff file and the wgsmaster.gbff'

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


class StrainParser(object):
    """Extract genes in nucleotide space."""

    def __init__(self):
        pass

    def run(self, genome_dir_file, out_dir, header):
        outf = open(os.path.join(out_dir, header + "_summary_file.txt"), "w")
        outf.write(
            "genome_id\tOrganism name\tstrain (assembly report)\tstrain (genomic gbff)\tstrain (wgsmaster gbff)\tcombined_strains\ttype_material\n")
        pattern = re.compile('^\s+\/strain=".+')
        with open(genome_dir_file, 'r') as genomelistfile:
            for line in genomelistfile:
                line_split = line.strip().split('\t')

                genome_id = line_split[0]
                genome_path = line_split[1]
                genome_dir_id = os.path.basename(os.path.normpath(genome_path))

                species = ''
                strain = ''
                strain_gff = ''
                typemat = ''
                wstrain = ''
                gstrains = []
                gstrain = ""

                genbank_file = os.path.join(
                    line_split[1], genome_dir_id + '_assembly_report.txt')
                if os.path.exists(genbank_file):
                    with open(genbank_file, 'r') as gfile:
                        for gline in gfile:
                            if gline.startswith('# Organism name: '):
                                speline = gline.replace("# Organism name:", "")
                                species = re.sub(
                                    r'\([^)]*\)', '', speline).strip()
                            if gline.startswith('# Infraspecific name:'):
                                strline = gline.replace(
                                    "# Infraspecific name:", "")
                                strain = strline.replace(
                                    "strain=", "").strip()
                            # if gline.startswith("# Relation to type material: assembly from type material"):
                            #    typemat = "True"
                            if gline.startswith("# Relation to type material:"):
                                processed_line = gline.replace(
                                    "# Relation to type material:", '')
                                typemat = processed_line.strip()
                            if typemat != "" and strain != "" and species != "":
                                break

                wgsmaster_file = os.path.join(
                    line_split[1], genome_dir_id + '_wgsmaster.gbff')
                if os.path.exists(wgsmaster_file):
                    with open(wgsmaster_file, 'r') as wfile:
                        for wline in wfile:
                            if pattern.match(wline):
                                wstrain = wline.replace(
                                    "/strain=", '').replace('"', '').rstrip().lstrip()

                genomic_file = os.path.join(
                    line_split[1], genome_dir_id + '_genomic.gbff')
                if os.path.exists(genomic_file):
                    with open(genomic_file, 'r') as gfile:
                        for gline in gfile:
                            if pattern.match(gline):
                                gstrain = gline.replace(
                                    "/strain=", '').replace('"', '').rstrip().lstrip()
                                gstrains.append(gstrain)

                if 'substr.' in strain:
                    print "{1} strain {0}".format(strain, genome_id)
                    strain = strain.split("substr.")[0].rstrip()
                if 'substr.' in gstrain:
                    print "{1} gstrain {0}".format(gstrain, genome_id)
                    gstrain = gstrain.split("substr.")[0].rstrip()
                if 'substr.' in wstrain:
                    print "{1} wstrain {0}".format(wstrain, genome_id)
                    wstrain = wstrain.split("substr.")[0].rstrip()

                outf.write("{3}\t{0}\t{1}\t{2}\t{5}\t{6}\t{4}\n".format(species, strain, gstrain,
                                                                        genome_id, typemat, wstrain, ' / '.join([strain, ','.join(set(gstrains)), wstrain])))
        outf.close()


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--genome_dir_file',
                        help='file indicating path to all genome directories')
    parser.add_argument('--out_dir', help='outputdirectory')
    parser.add_argument(
        '--header', help='header that will be concatenated to the output files (gb,rs,gbk,....)')

    args = parser.parse_args()

    try:
        p = StrainParser()
        p.run(args.genome_dir_file, args.out_dir, args.header)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
