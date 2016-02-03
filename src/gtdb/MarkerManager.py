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
from psycopg2.extensions import AsIs

import Tools
from Exceptions import GenomeDatabaseError


class MarkerManager(object):
    """Manages markers."""

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

    def externalMarkerIdsToMarkerIds(self, external_ids):
        """Get database marker identifiers from external identifiers.

        Parameters
        ----------
        external_ids : list
            List of external marker ids.

        Returns
        -------
        list
            List of database marker ids.
        """

        try:
            map_databases_to_ids = {}

            for external_id in external_ids:
                try:
                    (database_prefix, database_specific_id) = external_id.split(
                        "_", 1)
                except ValueError:
                    raise GenomeDatabaseError(
                        "All marker ids must have the form <prefix>_<id>. Offending id: %s" % str(external_id))

                if database_prefix not in map_databases_to_ids:
                    map_databases_to_ids[database_prefix] = {}
                map_databases_to_ids[database_prefix][
                    database_specific_id] = external_id

            temp_table_name = Tools.generateTempTableName()

            if len(map_databases_to_ids.keys()):
                self.cur.execute("CREATE TEMP TABLE %s (prefix text)" %
                            (temp_table_name,))
                query = "INSERT INTO {0} (prefix) VALUES (%s)".format(
                    temp_table_name)
                self.cur.executemany(query, [(x,)
                                        for x in map_databases_to_ids.keys()])
            else:
                raise GenomeDatabaseError(
                    "No marker databases found for these ids. %s" % str(external_ids))

            # Find any given database prefixes that arent in the marker
            # databases
            query = ("SELECT prefix FROM {0} " +
                     "WHERE prefix NOT IN ( " +
                     "SELECT external_id_prefix " +
                     "FROM marker_databases)").format(temp_table_name)

            self.cur.execute(query)

            missing_marker_sources = {}
            for (query_prefix,) in self.cur:
                missing_marker_sources[query_prefix] = map_databases_to_ids[
                    query_prefix].values()

            if len(missing_marker_sources.keys()):
                errors = []
                for (source_prefix, offending_ids) in missing_marker_sources.items():
                    errors.append("(%s) %s" %
                                  (source_prefix, str(offending_ids)))
                raise GenomeDatabaseError("Cannot find the relevant marker database id for the following ids, check the IDs are correct: " +
                                          ", ".join(errors))

            # All genome sources should be good, find ids
            result_ids = []
            for database_prefix in map_databases_to_ids.keys():

                # Create a table of requested external ids from this genome
                # source
                temp_table_name = Tools.generateTempTableName()
                self.cur.execute(
                    "CREATE TEMP TABLE %s (id_in_database text)" % (temp_table_name,))
                query = "INSERT INTO {0} (id_in_database) VALUES (%s)".format(
                    temp_table_name)
                self.cur.executemany(
                    query, [(x,) for x in map_databases_to_ids[database_prefix].keys()])

                # Check to see if there are any that don't exist
                query = ("SELECT id_in_database FROM {0} " +
                         "WHERE id_in_database NOT IN ( " +
                         "SELECT id_in_database " +
                         "FROM markers, marker_databases " +
                         "WHERE marker_database_id = marker_databases.id " +
                         "AND external_id_prefix = %s)").format(temp_table_name)

                self.cur.execute(query, (database_prefix,))

                missing_ids = []
                for (id_in_database,) in self.cur:
                    missing_ids.append(database_prefix + "_" + id_in_database)

                if missing_ids:
                    raise GenomeDatabaseError(
                        "Cannot find the the following marker ids, check the IDs are correct: %s" % str(missing_ids))

                # All exist, so get their ids.
                query = ("SELECT markers.id FROM markers, marker_databases " +
                         "WHERE marker_database_id = marker_databases.id " +
                         "AND id_in_database IN ( " +
                         "SELECT id_in_database " +
                         "FROM {0} )" +
                         "AND external_id_prefix = %s").format(temp_table_name)

                self.cur.execute(query, (database_prefix,))

                for (marker_id,) in self.cur:
                    result_ids.append(marker_id)

        except GenomeDatabaseError as e:
            raise e

        return result_ids

    def printMarkerDetails(self, marker_id_list):
        """Print marker gene details.

        Parameters
        ----------
        marker_id_list : iterable
            Unique identifier of markers in database.

        Returns
        -------
        bool
            True if successful, else False.
        """

        try:
            if not marker_id_list:
                raise GenomeDatabaseError(
                    "Unable to print markers. No markers found.")

            columns = "markers.id, markers.name, description, " + \
                "external_id_prefix || '_' || id_in_database as external_id, size"

            self.cur.execute("SELECT " + columns + " FROM markers " +
                                "LEFT OUTER JOIN users ON markers.owner_id = users.id " +
                                "JOIN marker_databases AS databases ON marker_database_id = databases.id " +
                                "AND markers.id in %s " +
                                "ORDER BY markers.id ASC", (tuple(marker_id_list),))

            # print table
            header = ("Marker ID", "Name", "Description", "Length (aa)")

            rows = []
            for (_marker_id, name, description, external_id, size) in self.cur:
                rows.append((external_id, name, description, size))

        except GenomeDatabaseError as e:
            raise e

        return header, rows
