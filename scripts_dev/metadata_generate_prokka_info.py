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

__prog_name__ = 'metadata_generate_prokka_info.py'
__prog_desc__ = ('Parse prokka_report.log file to generate metadata for the metadata_ncbi table.')

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2016'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@qfab.org'
__status__ = 'Development'


import os
import argparse
import sys


class ProkkaParser(object):

    def __init__(self):
        self.log = 1

    def generateMetadata(self, metadata_file, prokka_file, outfile):

        with open(metadata_file) as f:
            f.readline()
            genome_ids = [line.split(',')[0] for line in f]
        print genome_ids

        prokka_ids = {}
        with open(prokka_file) as f:
            for line in f:
                raw_prokka_id = line.split('\t')[1]
                if raw_prokka_id.startswith('GCA_'):
                    prokka_id = 'GB_' + raw_prokka_id
                elif raw_prokka_id.startswith('GCF_'):
                    prokka_id = 'RS_' + raw_prokka_id
                else:
                    sys.exit('Header is unknown.')
                status = line.split('\t')[2]
                prokka_ids[prokka_id] = status

        file_writer = open(outfile, 'w')
        file_writer.write("genome_id\tncbi_called_genes\n")
        for k, v in prokka_ids.iteritems():
            if k in genome_ids:
                file_writer.write("{0}\t{1}\n".format(k, 'TRUE'))
        file_writer.close()


if __name__ == "__main__":
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--prokka_report', dest="prokka_file", required=True, help='Report listing Genomes where genes have been called using Prokka')
    parser.add_argument('--metadata_file', dest="meta_file", required=True, help='Metadata files coming from GTDB metadata export')
    parser.add_argument('--output_file', dest="outfile", required=True, help='Output file')

    args = parser.parse_args()

    try:
        prokka_parser = ProkkaParser()
        prokka_parser.generateMetadata(args.meta_file, args.prokka_file, args.outfile)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
