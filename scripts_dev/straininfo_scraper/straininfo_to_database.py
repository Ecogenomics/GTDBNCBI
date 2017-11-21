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

__prog_name__ = 'update_database_from_ftp.py'
__prog_desc__ = ('Update the STRAININFO tables in GTDB. ' +
                 'STRAININFO tables are independant of the metadata GTDB information')

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

from database_configuration import GenomeDatabaseConnectionSTINFOUpdate


class UpdateSTRAININFODatabase(object):

    def __init__(self, path):

        self.path = path
        self.straininfo_strains_file = os.path.join(path, 'straininfo_strains.tsv')

        self.temp_con = GenomeDatabaseConnectionSTINFOUpdate.GenomeDatabaseConnectionSTINFOUpdate()
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

    def runUpdate(self):
        # Check if the files exist:
        if os.path.isfile(self.straininfo_strains_file):
            self.temp_cur.execute('TRUNCATE straininfo_strains;')
            print "Deletion straininfo_strains done"
            fr = open(self.straininfo_strains_file)
            fr.readline()
            self.temp_cur.copy_from(fr, 'straininfo_strains')
            print 'Copy straininfo_strains done'
            self.temp_con.commit()

        else:
            print 'Some files are missing in {0}'.format(self.path)
        self.temp_con.ClosePostgresConnection()


if __name__ == "__main__":
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--straininfo_dir', dest="straininfo_dir",
                        required=True, help='Directory to the STRAININFO files straininfo_strains.tsv')

    args = parser.parse_args()

    try:
        update_straininfo_mngr = UpdateSTRAININFODatabase(args.straininfo_dir)
        update_straininfo_mngr.runUpdate()
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
