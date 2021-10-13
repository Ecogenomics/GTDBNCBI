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

__prog_name__ = 'metadata_add_representatives_to_database.py'
__prog_desc__ = 'Add representative genomes to database.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2019'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import csv
import logging
import argparse
import tempfile
from collections import defaultdict

from biolib.common import canonical_gid

from gtdb import GenomeDatabase
from gtdb.Exceptions import (GenomeDatabaseError, 
                                DumpDBErrors, 
                                DumpDBWarnings, 
                                ErrorReport)
                                
import psycopg2
from psycopg2.extensions import AsIs


class AddRepresentativeGenomes(object):
  """Populate 'gtdb_genome_representative' and 'gtdb_representative' fields in database."""

  def __init__(self):    
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

  def run(self, final_cluster_file, gtdb_version):
    """Add metadata."""

    # clear representative fields
    self.logger.info('Setting GTDB representative fields to NULL.')
    self.setup_db(gtdb_version)
    cur = self.db.conn.cursor()
        
    q = ("UPDATE metadata_taxonomy SET gtdb_representative = NULL, gtdb_genome_representative = NULL")
    cur.execute(q)
    self.db.conn.commit()
    
    # mark all genomes as not being representatives and get translation
    # between canonical genome IDs and NCBI accessions
    q = ("SELECT accession FROM metadata_view")
    cur.execute(q)
    
    gid_to_ncbi_accn = {}
    is_rep = {}
    for r in cur:
        ncbi_accn = r[0]
        gid_to_ncbi_accn[canonical_gid(ncbi_accn)] = ncbi_accn
        is_rep[ncbi_accn] = False
    
    # determine representative assignment of genomes
    temp_genome_rep_file = tempfile.NamedTemporaryFile(delete=False)
    num_sp_reps = 0
    with open(final_cluster_file) as f:
        headers = f.readline().strip().split('\t')
        
        rep_index = headers.index('Representative')
        clustered_genomes_index = headers.index('Clustered genomes')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            rep_accn = gid_to_ncbi_accn[line_split[rep_index]]
            if len(line_split) > clustered_genomes_index:
                gids = [gid.strip() for gid in line_split[clustered_genomes_index].split(',')]
                for gid in gids:
                    ncbi_accn = gid_to_ncbi_accn[gid]
                    temp_genome_rep_file.write('%s\t%s\n' % (
                                                    ncbi_accn, 
                                                    rep_accn))
                
            temp_genome_rep_file.write('%s\t%s\n' % (
                                        rep_accn, 
                                        rep_accn))
            is_rep[rep_accn] = True
            num_sp_reps += 1
    temp_genome_rep_file.close()
    
    print('Identified %d species clusters.' % num_sp_reps)
    print('Identified %d genomes marked as representatives.' % sum([1 for rid in is_rep if is_rep[rid]]))
     
    cmd = 'gtdb -r metadata import --table metadata_taxonomy --field gtdb_genome_representative --type TEXT --metadatafile %s' % (temp_genome_rep_file.name)
    print(cmd)
    os.system(cmd)
    os.remove(temp_genome_rep_file.name)
    
    # mark representative genomes
    temp_rep_file = tempfile.NamedTemporaryFile(delete=False)
    for rep_accn, rep_status in list(is_rep.items()):
        temp_rep_file.write('%s\t%s\n' % (rep_accn, str(rep_status)))
    temp_rep_file.close()
    
    cmd = 'gtdb -r metadata import --table metadata_taxonomy --field gtdb_representative --type BOOLEAN --metadatafile %s' % (temp_rep_file.name)
    print(cmd)
    os.system(cmd)
    os.remove(temp_rep_file.name)  
    cur.close()

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('final_cluster_file', help="clusters for named species")
    parser.add_argument('gtdb_version', help='GTDB database version (i.e., gtdb_releaseX)')

    args = parser.parse_args()

    try:
        p = AddRepresentativeGenomes()
        p.run(args.final_cluster_file, args.gtdb_version)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
