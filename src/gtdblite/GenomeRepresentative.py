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

import logging

import psycopg2


class GenomeRepresentative(object):
    ''''Manage genome representatives.'''

    def __init__(self, cur, currentUser):
        """Initialize.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        currentUser : User
            Current user of database.
        """

        self.logger = logging.getLogger()

        self.cur = cur
        self.currentUser = currentUser

    def determineRepresentative(self):
        pass

    def __assignRepresentatives(self):
        """Assign genomes to a representative.

        This function will try to assign all
        non-representative genomes to a
        representative. In normal use, this
        function is not called, but it is helpful
        when large numbers of genomes are added
        directly to the database.
        """

        # get representative genomes
        self.cur.execute("SELECT id "
                         "FROM genomes "
                         "WHERE gtdb_representative = 'TRUE'")
        rep_genome_ids = [genome_id for genome_id in self.cur.fetchall()]
        print 'rep_genome_ids', len(rep_genome_ids)

        # compare genomes to representatives
        self.cur.execute("SELECT id FROM genomes")
        for genome_id in self.cur:
            for rep_genome_id in rep_genome_ids:
                self.cur.execute("SELECT id "
                         "FROM genomes "
                         "WHERE gtdb_representative = 'TRUE'")

