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


__prog_name__ = 'silva_aggregate_blast_hits.py'
__prog_desc__ = 'Aggregates all hits from all blast files into a single file.'

__author__ = 'Aaron Mussig'
__copyright__ = 'Copyright 2018'
__credits__ = ['Aaron Mussig']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Aaron Mussig'
__email__ = 'a.mussig@uq.edu.au'
__status__ = 'Development'

import argparse
import os
import sys
from multiprocessing import Pool

import Config
import GenomeDatabase
from Exceptions import GenomeDatabaseError, ErrorReport, DumpDBErrors


def aggregate_worker(worker_args):
    """A worker for multi-processing."""
    gid, full_path = worker_args
    output = str()
    if os.path.isfile(full_path):
        with open(full_path, 'r') as f:
            f.readline()
            for line in f.readlines():
                output += str(gid) + '\t' + line
    return output


class Aggregate(object):

    def __int__(self):
        pass

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

    def run(self, gtdb_release, silva_file_name, threads, output_file):
        # type: (str, str, int, str) -> None
        """Visit each genome directory in the specified release and extract the SILVA summary file."""

        self.setup_db(gtdb_release)
        cur = self.db.conn.cursor()

        cur.execute("""SELECT mt.id, gs.name, g.fasta_file_location FROM metadata_taxonomy mt 
                       INNER JOIN genomes g on mt.id = g.id INNER JOIN genome_sources gs on g.genome_source_id = gs.id
                       WHERE gtdb_representative""")

        queue = list()
        for gid, source, path in cur.fetchall():
            if source == 'user':
                root_path = Config.GTDB_GENOME_USR_DIR
            elif source == 'GenBank':
                root_path = Config.GTDB_GENOME_GBK_DIR
            elif source == 'RefSeq':
                root_path = Config.GTDB_GENOME_RSQ_DIR
            full_path = root_path + '/'.join(path.split('/')[0:-1]) + '/rna_silva/' + silva_file_name
            queue.append((gid, full_path))

        p = Pool(threads)
        results = p.map(aggregate_worker, queue)

        with open(output_file, 'w') as f:
            f.write(''.join(results))


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_release',
                        help='The release version of GTDB (e.g. gtdb_release86).')
    parser.add_argument('silva_file_name',
                        help='The format of the SILVA summary file (ssu.taxonomy.tsv, lsu_23S.taxonomy.tsv).')
    parser.add_argument('threads', type=int, help='The number of threads to use.')
    parser.add_argument('output_file', help='The output file to write the aggregated data to (TSV).')
    args = parser.parse_args()

    try:
        p = Aggregate()
        p.run(args.gtdb_release, args.silva_file_name, args.threads, args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
