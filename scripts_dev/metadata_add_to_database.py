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

__prog_name__ = 'metadata_add_to_database.py'
__prog_desc__ = 'Add columns in a metadata file to the GTDB.'

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
import logging
import argparse
import tempfile
from collections import defaultdict

from gtdb import GenomeDatabase
from gtdb.Exceptions import (GenomeDatabaseError,
                             DumpDBErrors,
                             DumpDBWarnings,
                             ErrorReport)

import psycopg2
from psycopg2.extensions import AsIs


class AddMetadata(object):
    """Add all columns in a metadata file to the GTDB."""

    def __init__(self):
        """Initializatin."""
        logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S",
                            level=logging.DEBUG)
        self.logger = logging

    def setup_db(self, gtdb_version):
        """Setup database."""

        # initialise the backend
        self.db = GenomeDatabase.GenomeDatabase(1, False)
        self.db.conn.MakePostgresConnection(gtdb_version)

        # login
        try:
            self.db.Login(None, False)
        except GenomeDatabaseError as e:
            self.db.conn.ClosePostgresConnection()
            ErrorReport(e.message + " The following error(s) were reported:\n")
            DumpDBErrors(self.db)
            sys.exit(-1)

    def run(self, metadata_file,
            metadata_desc_file,
            genome_list_file,
            do_not_null_field,
            gtdb_version):
        """Add metadata."""

        # get fields in metadata file
        with open(metadata_file) as f:
            metadata_fields = f.readline().strip().split('\t')[1:]
        self.logger.info(
            'Metadata file contains {} fields.'.format(len(metadata_fields)))
        self.logger.info('Fields: %s' % ', '.join(metadata_fields))

        # get database table and data type of each metadata field
        metadata_type = {}
        metadata_table = {}
        with open(metadata_desc_file) as f:
            for line in f:
                line_split = line.strip('\n').split('\t')
                field = line_split[0]
                if field in metadata_fields:
                    metadata_type[field] = line_split[2]
                    metadata_table[field] = line_split[3]
        self.logger.info('Identified {} matching fields in metadata description file.'.format(
            len(metadata_table)))
        self.logger.info('Fields: %s' % ', '.join(metadata_table))

        # set fields to NULL if requested
        if not do_not_null_field:
            response = ''
            while response.lower() not in ['y', 'n']:
                response = raw_input(
                    "Set fields to NULL for all genomes [y/n]: ")

            if response.lower() == 'y':
                self.logger.info('Connecting to %s.' % gtdb_version)
                self.setup_db(gtdb_version)
                cur = self.db.conn.cursor()
        
                for field in metadata_table:
                    q = ("UPDATE {} SET {} = NULL".format(
                        metadata_table[field], field))
                    print(q)
                    cur.execute(q)
                self.db.conn.commit()

                cur.close()
            elif response.lower() == 'n':
                pass
            else:
                self.logger.error('Unrecognized input.')
                sys.exit(-1)

        # get genomes to process
        genome_list = set()
        if genome_list_file:
            for line in open(genome_list_file):
                if len(line.split('\t')) >= len(line.split(',')):
                    genome_list.add(line.rstrip().split('\t')[0])
                else:
                    genome_list.add(line.rstrip().split(',')[0])

        # read metadata file
        metadata = defaultdict(lambda: defaultdict(str))
        with open(metadata_file) as f:
            fields = [x.strip() for x in f.readline().split('\t')]

            for line in f:
                line_split = line.rstrip('\n').split('\t')

                genome_id = line_split[0]
                # print line_split
                for i, value in enumerate(line_split[1:]):
                    metadata[fields[i + 1]][genome_id] = value

        # add each field to the database
        for field in metadata:
            temp_file = tempfile.NamedTemporaryFile(delete=False)

            if field not in metadata_type:
                continue

            data_type = metadata_type[field]
            table = metadata_table[field]

            records_to_update = 0
            for orig_genome_id, value in metadata[field].iteritems():

                try:
                    if float(value) and data_type in ['INT', 'INTEGER']:
                        # assume specified data type is correct and that we may need
                        # to cast floats to integers
                        value = str(int(float(value)))
                except:
                    pass

                if value.strip():
                    genome_id = str(orig_genome_id)
                    if genome_id.startswith('GCA_'):
                        genome_id = 'GB_' + genome_id
                    elif genome_id.startswith('GCF_'):
                        genome_id = 'RS_' + genome_id

                    if (not genome_list 
                            or genome_id in genome_list 
                            or orig_genome_id in genome_list):
                        temp_file.write('%s\t%s\n' % (genome_id, value))
                        records_to_update += 1

            temp_file.close()

            self.logger.info('Updating {} for {} genomes.'.format(
                field, records_to_update))

            cmd = 'gtdb -r metadata import --table %s --field %s --type %s --metadatafile %s' % (
                table, field, data_type, temp_file.name)
            print(cmd)
            os.system(cmd)
            os.remove(temp_file.name)


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'metadata_file', help='tab-separated values table with metadata')
    parser.add_argument(
        'metadata_desc_file', help='tab-separated values table with description of metadata fields')
    parser.add_argument(
        '--genome_list', help='only process genomes in this list', default=None)
    parser.add_argument('--do_not_null_field',
                        help='Do not set fields to NULL for all genomes before updating values',
                        action='store_true')
    parser.add_argument(
        '--gtdb_version', help='GTDB database version (i.e., gtdb_releaseX)')

    args = parser.parse_args()

    try:
        p = AddMetadata()
        p.run(args.metadata_file,
              args.metadata_desc_file,
              args.genome_list,
              args.do_not_null_field,
              args.gtdb_version)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
