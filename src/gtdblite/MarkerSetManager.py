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

        self.bacCanonicalMarkerSetId = 1
        self.arCanonicalMarkerSetId = 2

    def _confirm(self, msg):
        raw = raw_input(msg + " (y/N): ")
        if raw.upper() == "Y":
            return True
        return False

    def canonicalBacterialMarkers(self):
        """Get identifiers of canonical bacterial markers."""

        return self.getMarkerIdsFromMarkerSetIds([self.bacCanonicalMarkerSetId])

    def canonicalArchaealMarkers(self):
        """Get identifiers of canonical archaeal markers."""

        return self.getMarkerIdsFromMarkerSetIds([self.arCanonicalMarkerSetId])

    def concatenatedAlignedMarkers(self, db_genome_id, marker_id_index):
        """Create concatenated alignment of markers for genome.

        Markers are concatenated in the order specified by
        marker_id_index.

        Parameters
        ----------
        db_genome_id : str
            Unique database identifier of genome of interest.
        marker_id_index : d[marker_id] -> index
            Dictionary indicating the order to concatenate markers.

        Returns
        -------
        str
            Concatenated alignment.
        """

        query = ("SELECT marker_id, sequence " +
                 "FROM aligned_markers " +
                 "WHERE genome_id = %s AND marker_id = ANY(%s)")
        self.cur.execute(query, (db_genome_id, marker_id_index.keys()))

        concatenated_align = [None] * len(marker_id_index)
        for marker_id, sequence in self.cur:
            index = marker_id_index[marker_id]
            concatenated_align[index] = sequence

        return ''.join(concatenated_align)

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
            self.cur.executemany(
                query, [(marker_set_id, x) for x in marker_id_list])

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

            update_query = []
            params = []

            if name is not None:
                update_query.append("name = %s")
                params.append(name)

            if description is not None:
                update_query.append("description = %s")
                params.append(description)

            if private is not None:
                update_query.append("private = %s")
                params.append(private)

            if params:
                self.cur.execute("UPDATE marker_sets SET " + ",".join(update_query) +
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

                    query_is_empty = ("SELECT count(msc.marker_id) from marker_sets as ms " +
                                      "LEFT JOIN  marker_set_contents as msc on msc.set_id = ms.id " +
                                      "WHERE ms.id = {0} " +
                                      "GROUP BY ms.id").format(marker_set_id)
                    self.cur.execute(query_is_empty)

                    count = self.cur.fetchone()
                    if count[0] == 0:
                        # We delete the list because it's empty
                        query_del_set = ("DELETE FROM marker_sets WHERE id = {0} ").format(
                            marker_set_id)
                        self.cur.execute(query_del_set)

                else:
                    raise GenomeDatabaseError(
                        "Unknown marker set edit operation: %s" % operation)

        except GenomeDatabaseError as e:
            raise e

        return True

    def deleteMarkerSets(self, marker_set_ids):
        """Delete marker set and associated markers.

        Parameters
        ----------
        marker_set_ids : iterable
            Unique identifier of marker sets in database.

        Returns
        -------
        bool
            True if successful.
        """

        for marker_set_id in marker_set_ids:
            try:
                edit_permission = self.permissionToModify(marker_set_id)
                if edit_permission is False:
                    raise GenomeDatabaseError(
                        "Insufficient permissions to delete marker set. Offending marker set id: {0}".format(marker_set_id))

                if not self._confirm("Are you sure you want to delete {0} set(s) (this action cannot be undone)".format(len(marker_set_ids))):
                    raise GenomeDatabaseError("User aborted database action.")

                list_marker_ids = self.getMarkerIdsFromMarkerSetIds(
                    [marker_set_id])
                self.editMarkerSet(marker_set_id, list_marker_ids, 'remove')
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

    def getMarkerIdsFromMarkerSetIds(self, marker_set_ids):
        """Get marker identifiers within specific marker set.

        Parameters
        ----------
        marker_set_ids : iterable
            Identifiers of marker sets.

        Returns
        -------
        list
            Identifier of markers in marker sets.
        """

        self.cur.execute("SELECT id, owner_id, owned_by_root, private " +
                         "FROM marker_sets " +
                         "WHERE id = %s ", (tuple(marker_set_ids),))

        result = self.cur.fetchone()

        if not result:
            self.ReportError(
                "At least one marker set is invalid: %s" % str(marker_set_ids))
            return None

        self.cur.execute("SELECT marker_id " +
                         "FROM marker_set_contents " +
                         "WHERE set_id = %s ", (tuple(marker_set_ids),))

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
