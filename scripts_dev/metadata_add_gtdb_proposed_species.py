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

__prog_name__ = 'metadata_add_gtdb_proposed_species.py'
__prog_desc__ = 'Add proposed GTDB species to database.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2019'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
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

from gtdb import GenomeDatabase
from gtdb.Exceptions import (GenomeDatabaseError, 
                                DumpDBErrors, 
                                DumpDBWarnings, 
                                ErrorReport)
                                
import psycopg2
from psycopg2.extensions import AsIs


class AddSpecies(object):
  """Populate 'gtdb_species' field in database."""

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

  def run(self, user_cluster_file, gtdb_version):
    """Add metadata."""
    
    self.logger.info('Connecting to %s.' % gtdb_version)
    self.setup_db(gtdb_version)
    
    # clear GTDB species field
    cur = self.db.conn.cursor()
        
    q = ("UPDATE metadata_taxonomy SET gtdb_proposed_species = NULL")
    cur.execute(q)
    self.db.conn.commit()

    # set proposed species
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    with open(user_cluster_file) as f:
        header = f.readline().strip().split('\t')
        
        gid_index = header.index('Type genome')
        sp_index = header.index('NCBI species')
        cluster_index = header.index('Clustered genomes')
        
        for line in f:
            line_split = line.strip('\n').split('\t')
            
            sp = line_split[sp_index]
            
            rep_id = line_split[gid_index]
            temp_file.write('%s\t%s\n' % (rep_id, sp))
            
            if line_split[cluster_index]:
                for gid in line_split[cluster_index].split(','):
                    temp_file.write('%s\t%s\n' % (gid, sp))
        
    temp_file.close()
    
    cmd = 'gtdb -r metadata import --table metadata_taxonomy --field gtdb_proposed_species --type TEXT --metadatafile %s' % (temp_file.name)
    print cmd
    os.system(cmd)
    os.remove(temp_file.name)  
    cur.close()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('user_cluster_file', help="GTDB species clusters for all NCBI and User genomes")
    parser.add_argument('gtdb_version', help='GTDB database version (i.e., gtdb_releaseX)')

    args = parser.parse_args()

    try:
        p = AddSpecies()
        p.run(args.user_cluster_file, args.gtdb_version)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
