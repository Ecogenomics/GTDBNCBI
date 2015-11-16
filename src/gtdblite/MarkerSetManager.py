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


class MarkerSetManager(object):
    """Manages marker sets."""

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

    def createMarkerSet(self, marker_id_list, name, description, owner_id=None, private=None):
        """Create a new marker set.

        Parameters
        ----------
        marker_id_list : list
            A list of marker ids to add to the new list.
        name : str
            Name of the newly created list.
        description : str
            Description of the newly created list.
        owner_id : str
            The id of the user who will own this list.
        private : bool
            Denotes whether this list is public or private.

        Returns
        -------
        int
            Database identifier of newly created marker list.
        """
        try:
            if (owner_id is None):
                if not self.currentUser.isRootUser():
                    raise GenomeDatabaseError(
                        "Only the root user can create root owned lists.")
            else:
                if (not self.currentUser.isRootUser()) and (self.currentUser.getUserId() != owner_id):
                    raise GenomeDatabaseError(
                        "Only the root user may create sets on behalf of other people.")

            query = "INSERT INTO marker_sets (name, description, owned_by_root, owner_id, private) VALUES (%s, %s, %s, %s, %s) RETURNING id"
            self.cur.execute(
                query, (name, description, owner_id is None, owner_id, private))
            (marker_set_id,) = self.cur.fetchone()

            query = "INSERT INTO marker_set_contents (set_id, marker_id) VALUES (%s, %s)"
            self.cur.executemany(query, [(marker_set_id, x) for x in marker_id_list])

        except GenomeDatabaseError as e:
            raise e

        return marker_set_id

    def editMarkerSet(self, marker_set_id, marker_ids=None, operation=None, name=None, description=None, private=None):
        """Edit an existing marker set in the database.

        Parameters
        ----------
        marker_set_id : int
            Identifier of marker set in database.
        batchfile : str
            Filename of batch file describing markers to modify.
        marker_external_ids : list
            List of markers to modify.
        operation : str
            Operation to perform on genome list (add or remove).
        name : str
            Name of the newly created list.
        description : str
            Description of the newly created list.
        private : bool
            Denotes whether this list is public or private.

        Returns
        -------
        bool
            True if successful, else False
        """

        try:
            edit_permission = self.permissionToModify(marker_set_id)
            if edit_permission is None:
                raise GenomeDatabaseError(
                    "Unable to retrieve marker set id for editing. Offending set id: %s" % marker_set_id)
            elif edit_permission is False:
                raise GenomeDatabaseError(
                    "Insufficient permissions to edit this marker set. Offending set id: %s" % marker_set_id)

            update_query = ""
            params = []

            if name is not None:
                update_query += "name = %s"
                params.append(name)

            if description is not None:
                update_query += "description = %s"
                params.append(description)

            if private is not None:
                update_query += "private = %s"
                params.append(private)

            if params:
                self.cur.execute("UPDATE marker_sets SET " + update_query +
                            " WHERE id = %s", params + [marker_set_id])

            temp_table_name = Tools.generateTempTableName()

            if operation is not None:

                if len(marker_ids) == 0:
                    raise GenomeDatabaseError(
                        "No marker ids given to perform '%s' operation." % operation)

                self.cur.execute("CREATE TEMP TABLE %s (id integer)" %
                            (temp_table_name,))
                query = "INSERT INTO {0} (id) VALUES (%s)".format(
                    temp_table_name)
                self.cur.executemany(query, [(x,) for x in marker_ids])

                if operation == 'add':
                    query = ("INSERT INTO marker_set_contents (set_id, marker_id) " +
                             "SELECT %s, id FROM {0} " +
                             "WHERE id NOT IN ( " +
                             "SELECT marker_id " +
                             "FROM marker_set_contents " +
                             "WHERE set_id = %s)").format(temp_table_name)
                    self.cur.execute(query, (marker_set_id, marker_set_id))
                elif operation == 'remove':
                    query = ("DELETE FROM marker_set_contents " +
                             "WHERE set_id = %s " +
                             "AND marker_id IN ( " +
                             "SELECT id " +
                             "FROM {0})").format(temp_table_name)
                    self.cur.execute(query, [marker_set_id])
                else:
                    raise GenomeDatabaseError(
                        "Unknown marker set edit operation: %s" % operation)

        except GenomeDatabaseError as e:
            raise e

        return True

    def getAllMarkerIds(self):
        """Get identifiers for all markers.

        Returns
        -------
        list
            Database identifiers for all markers.
        """

        try:
            query = "SELECT id FROM markers"
            self.cur.execute(query)

            result_ids = []
            for (marker_id,) in self.cur:
                result_ids.append(marker_id)

        except GenomeDatabaseError as e:
            raise e

        return result_ids

    def getMarkerIdListFromMarkerSetId(self, marker_set_id):
        """Get marker identifiers within specific marker set.

        Parameters
        ----------
        marker_set_id : int
            Identifier of marker set.

        Returns
        -------
        list
            Identifier of markers in marker set.
        """

        self.cur.execute("SELECT id, owner_id, owned_by_root, private " +
                    "FROM marker_sets " +
                    "WHERE id = %s ", (tuple(marker_set_id),))

        result = self.cur.fetchone()

        if not result:
            self.ReportError("No marker set with id: %s" % str(marker_set_id))
            return None

        self.cur.execute("SELECT marker_id " +
                    "FROM marker_set_contents " +
                    "WHERE set_id = %s ", (tuple(marker_set_id),))

        return [marker_id for (marker_id,) in self.cur.fetchall()]

    def printMarkerSetsDetails(self, marker_set_ids):
        """Print marker set details.

        Parameters
        ----------
        genome_list_ids : iterable
            Unique identifier of marker sets in database.

        Returns
        -------
        list
            Column headers.
        list
            Content for each row.
        """

        try:
            if not marker_set_ids:
                raise GenomeDatabaseError(
                    "Unable to print marker set details: No marker sets given.")

            if not self.currentUser.isRootUser():
                self.cur.execute("SELECT id " +
                            "FROM marker_sets as sets " +
                            "WHERE sets.private = True " +
                            "AND sets.id in %s ", (tuple(marker_set_ids),))

                unviewable_set_ids = [set_id for (set_id,) in self.cur]
                if unviewable_set_ids:
                    raise GenomeDatabaseError(
                        "Insufficient privileges to view marker sets: %s." % str(unviewable_set_ids))

            self.cur.execute(
                "SELECT sets.id, sets.name, sets.description, sets.owned_by_root, users.username, count(contents.set_id) " +
                "FROM marker_sets as sets " +
                "LEFT OUTER JOIN users ON sets.owner_id = users.id " +
                "JOIN marker_set_contents as contents ON contents.set_id = sets.id " +
                "WHERE sets.id in %s " +
                "GROUP by sets.id, users.username " +
                "ORDER by sets.id asc ", (tuple(marker_set_ids),)
            )

            header = ("set_id", "name", "description", "owner", "marker_count")

            rows = []
            for (set_id, name, description, owned_by_root, username, marker_count) in self.cur:
                rows.append(
                    (set_id, name, description, ("root" if owned_by_root else username), marker_count))

        except GenomeDatabaseError as e:
            raise e

        return header, rows

    def permissionToModify(self, marker_set_id):
        """Check if user has permission to modify marker set.

        Parameters
        ----------
        marker_set_id : int
            Unique identifier of marker set in database.

        Returns
        -------
        bool
            True if has permission, else False.
        """

        try:
            self.cur.execute("SELECT owner_id, owned_by_root " +
                        "FROM marker_sets " +
                        "WHERE id = %s ", (marker_set_id,))

            result = self.cur.fetchone()

            if not result:
                raise GenomeDatabaseError(
                    "No marker set with id: %s" % str(marker_set_id))

            (owner_id, owned_by_root) = result

            if not self.currentUser.isRootUser():
                if owned_by_root or owner_id != self.currentUser.getUserId():
                    return False
            else:
                if not owned_by_root:
                    raise GenomeDatabaseError(
                        "Root user editing of other users marker sets not yet implmented.")

        except GenomeDatabaseError as e:
            raise e

        return True