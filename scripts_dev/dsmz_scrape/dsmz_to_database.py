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
__prog_desc__ = ('Update the DSMZ tables in GTDB. ' +
                 'DSMZ tables are independant of the metadata GTDB information')

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

from database_configuration import GenomeDatabaseConnectionDSMZUpdate


class UpdateDSMZDatabase(object):

    def __init__(self, path):

        self.path = path
        self.dsmz_genera_file = os.path.join(path, 'dsmz_genera.tsv')
        self.dsmz_strains_file = os.path.join(path, 'dsmz_strains.tsv')
        self.dsmz_species_file = os.path.join(path, 'dsmz_species.tsv')

        self.temp_con = GenomeDatabaseConnectionDSMZUpdate.GenomeDatabaseConnectionDSMZUpdate()
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

    def runUpdate(self):
        # Check if the files exist:
        if os.path.isfile(self.dsmz_genera_file) and os.path.isfile(self.dsmz_strains_file) and os.path.isfile(self.dsmz_species_file):
            self.temp_cur.execute('TRUNCATE dsmz_genera;')
            print "Deletion dsmz_genera done"
            fr = open(self.dsmz_genera_file)
            fr.readline()
            self.temp_cur.copy_from(fr, 'dsmz_genera')
            print 'Copy dsmz_genera done'
            self.temp_con.commit()

            self.temp_cur.execute('TRUNCATE dsmz_species;')
            print "Deletion dsmz_species done"
            fr = open(self.dsmz_species_file)
            fr.readline()
            self.temp_cur.copy_from(fr, 'dsmz_species')
            print 'Copy dsmz_species done'
            self.temp_con.commit()

            fr = open(self.dsmz_strains_file)
            fr.readline()
            self.temp_cur.execute('TRUNCATE dsmz_strains;')
            print "Deletion dsmz_strains done"

            self.temp_cur.copy_from(fr, 'dsmz_strains')
            print 'Copy dsmz_strains done'
            self.temp_con.commit()

        else:
            print 'Some files are missing in {0}'.format(self.path)
        self.temp_con.ClosePostgresConnection()


if __name__ == "__main__":
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--dsmz_dir', dest="dsmz_dir",
                        required=True, help='Directory to the DSMZ files dsmz_genera.tsv  dsmz_species.tsv  dsmz_strains.tsv')

    args = parser.parse_args()

    try:
        update_dsmz_mngr = UpdateDSMZDatabase(args.dsmz_dir)
        update_dsmz_mngr.runUpdate()
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
