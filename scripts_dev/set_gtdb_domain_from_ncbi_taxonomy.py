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

__prog_name__ = 'set_gtdb_domain_from_ncbi_domain.py'
__prog_desc__ = 'Set missing GTDB domain information to reflect NCBI domain.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2017'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import logging
from collections import defaultdict

from biolib.taxonomy import Taxonomy

from gtdb import GenomeDatabase
from gtdb.Exceptions import (GenomeDatabaseError, 
                                DumpDBErrors, 
                                DumpDBWarnings, 
                                ErrorReport)
                                
import psycopg2
from psycopg2.extensions import AsIs

class Script():
    def __init__(self):
        """Initialization."""
        
        logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S",
                            level=logging.DEBUG)
        self.logger = logging
        
    def setup_db(self, gtdb_version, threads):
        """Setup database."""

        # initialise the backend
        self.db = GenomeDatabase.GenomeDatabase(threads, False)
        self.db.conn.MakePostgresConnection(gtdb_version)

        # login
        try:
            self.db.Login(None, False)
        except GenomeDatabaseError as e:
            self.db.conn.ClosePostgresConnection()
            ErrorReport(e.message + " The following error(s) were reported:\n")
            DumpDBErrors(self.db)
            sys.exit(-1)
            
    def run(self):
        """Run script."""
        
        self.logger.info('Identifying NCBI genomes with missing domain information.')
        
        cur = self.db.conn.cursor()
        
        q = ("SELECT id, ncbi_taxonomy FROM metadata_taxonomy "
                + "WHERE (gtdb_domain IS NULL or gtdb_domain = 'd__') and ncbi_taxonomy IS NOT NULL")
        cur.execute(q)
        
        missing_domain_info = []
        for id, ncbi_taxonomy in cur:
            ncbi_domain = map(str.strip, ncbi_taxonomy.split(';'))[0]
            if ncbi_domain[0:3] != 'd__':
                self.logger.error('NCBI domain has the incorrect prefix: %s' % ncbi_domain)
                sys.exit()
             
            #gtdb_taxonomy = list(Taxonomy.rank_prefixes)
            #gtdb_taxonomy[0] = ncbi_domain
            #gtdb_taxonomy = ';'.join(gtdb_taxonomy)
            missing_domain_info.append([ncbi_domain, id])
 
        q = "UPDATE metadata_taxonomy SET gtdb_domain = %s WHERE id = %s"
        cur.executemany(q, missing_domain_info)
        
        self.db.conn.commit()
        cur.close()
                
        print 'NCBI genomes that were missing GTDB domain info: %d' % len(missing_domain_info)
 
if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_version', help='GTDB database version (i.e., gtdb_releaseX)')
    parser.add_argument('-t', '--threads', help='number of CPUs to use', type=int, default=32)

    args = parser.parse_args()

    try:
        p = Script()
        
        # setup database
        p.setup_db(args.gtdb_version, args.threads)
    
        p.run()
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise