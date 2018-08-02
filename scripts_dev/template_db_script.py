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

__prog_name__ = '<script name>.py'
__prog_desc__ = 'Template for scripts requiring direct access to GTDB database.'

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
        
    def setup_db(self, threads):
        """Setup database."""
        
        self.logger.info('Make sure you are running against the correct GTDB database!')
        
        # initialise the backend
        self.db = GenomeDatabase.GenomeDatabase(threads, False)
        self.db.conn.MakePostgresConnection()

        # login
        try:
            self.db.Login(None, False)
        except GenomeDatabaseError as e:
            self.db.conn.ClosePostgresConnection()
            ErrorReport(e.message + " The following error(s) were reported:\n")
            DumpDBErrors(self.db)
            sys.exit(-1)
            
    def run(self, argument, threads):
        """Run script."""
        
        cur = self.db.conn.cursor()
        
        q = ("SELECT id, gtdb_domain, ncbi_taxonomy FROM metadata_taxonomy")
        cur.execute(q)
        
        for id, gtdb_domain, ncbi_taxonomy in cur:
            print id, gtdb_domain, ncbi_taxonomy
            
        self.db.conn.commit()
        cur.close()
 
if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('argument', help="help message")
    parser.add_argument('-t', '--threads', help='number of CPUs to use', type=int, default=1)

    args = parser.parse_args() 

    try:
        p = Script()
        
        # setup database
        p._setup_db(args.threads)
    
        p.run(args.argument, args.threads)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise