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


class GenomeListManager(object):
    """Manages genome lists."""

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

    # TODO: This should not be here, technically the backend is agnostic so
    # shouldn't assume command line.
    def _confirm(self, msg):
        raw = raw_input(msg + " (y/N): ")
        if raw.upper() == "Y":
            return True
        return False

    def addGenomeList(self, genome_id_list, name, description, owner_id=None, private=None):
        """Creates a new genome list in the database.

        Parameters
        ----------
        genome_id_list : list
            A list of genome ids to add to the new list.
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
            Database identifier of newly created genome list.
        """

        try:
            if (owner_id is None):
                if not self.currentUser.isRootUser():
                    raise GenomeDatabaseError(
                        "Only the root user can create root owned lists.")
            else:
                if (not self.currentUser.isRootUser()) and (self.currentUser.getUserId() != owner_id):
                    raise GenomeDatabaseError(
                        "Only the root user may create lists on behalf of other people.")

            query = "INSERT INTO genome_lists (name, description, owned_by_root, owner_id, private) VALUES (%s, %s, %s, %s, %s) RETURNING id"
            self.cur.execute(
                query, (name, description, owner_id is None, owner_id, private))
            (genome_list_id,) = self.cur.fetchone()

            query = "INSERT INTO genome_list_contents (list_id, genome_id) VALUES (%s, %s)"
            self.cur.executemany(query, [(genome_list_id, x)
                                         for x in genome_id_list])
        except GenomeDatabaseError as e:
            raise e

        return genome_list_id

    def editGenomeList(self, genome_list_id, genome_ids=None, operation=None, name=None, description=None, private=None):
        """Edit an existing genome list in the database.

        Parameters
        ----------
        genome_list_id : int
            Identifier of genome list in database.
        genome_ids : list
            A list of genome ids to be modified.
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
            edit_permission = self.permissionToModify(genome_list_id)
            if edit_permission is None:
                raise GenomeDatabaseError(
                    "Unable to retrieve genome list id for editing. Offending list id: %s" % genome_list_id)
            if edit_permission is False:
                raise GenomeDatabaseError(
                    "Insufficient permissions to edit this genome list. Offending list id: %s" % genome_list_id)

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
                self.cur.execute("UPDATE genome_lists SET " + ",".join(update_query) +
                                 " WHERE id = %s", params + [genome_list_id])

            temp_table_name = Tools.generateTempTableName()

            if operation is not None:
                if len(genome_ids) == 0:
                    raise GenomeDatabaseError(
                        "No genome ids given to perform '%s' operation." % operation)

                self.cur.execute("CREATE TEMP TABLE %s (id integer)" %
                                 (temp_table_name,))
                query = "INSERT INTO {0} (id) VALUES (%s)".format(
                    temp_table_name)
                self.cur.executemany(query, [(x,) for x in genome_ids])

                if operation == 'add':
                    query = ("INSERT INTO genome_list_contents (list_id, genome_id) " +
                             "SELECT %s, id FROM {0} " +
                             "WHERE id NOT IN ( " +
                             "SELECT genome_id " +
                             "FROM genome_list_contents " +
                             "WHERE list_id = %s)").format(temp_table_name)
                    self.cur.execute(query, (genome_list_id, genome_list_id))
                elif operation == 'remove':
                    query = ("DELETE FROM genome_list_contents " +
                             "WHERE list_id = %s " +
                             "AND genome_id IN ( " +
                             "SELECT id " +
                             "FROM {0})").format(temp_table_name)
                    self.cur.execute(query, [genome_list_id])

                    query_is_empty = ("SELECT count(glc.genome_id) from genome_lists as gl " +
                                      "LEFT JOIN  genome_list_contents as glc on glc.list_id = gl.id " +
                                      "WHERE gl.id = {0} " +
                                      "GROUP BY gl.id").format(genome_list_id)
                    self.cur.execute(query_is_empty)

                    count = self.cur.fetchone()

                    if count[0] == 0:
                        # We delete the list because it's empty
                        query_del_list = ("DELETE FROM genome_lists WHERE id = {0} ").format(
                            genome_list_id)
                        self.cur.execute(query_del_list)

                else:
                    raise GenomeDatabaseError(
                        "Unknown genome set edit operation: %s" % operation)
        except GenomeDatabaseError as e:
            raise e

        return True

    def getGenomeIdsFromGenomeListIds(self, genome_list_ids):
        '''
        Function: getGenomeIdsFromGenomeListIds
        Given a list of ids, return all the ids of the genomes contained

        :param genome_list_ids: A list of genome list ids whose contents needs to be retrieved.

        Returns:
        A list of all the genome ids contained within the specified genome
        list(s), False on failure.
        '''
        try:

            temp_table_name = Tools.generateTempTableName()

            if genome_list_ids:
                try:
                    self.cur.execute("CREATE TEMP TABLE %s (id integer)" %
                                     (temp_table_name,))
                    query = "INSERT INTO {0} (id) VALUES (%s)".format(
                        temp_table_name)
                    self.cur.executemany(query, [(genome_list_id,) for genome_list_id in genome_list_ids])
                except:
                    raise GenomeDatabaseError("Unable to process genome list ID: %s" % str(genome_list_id))
            else:
                raise GenomeDatabaseError(
                    "No genome lists given. Cannot retrieve IDs")

            # Find any ids that don't have genome lists
            query = ("SELECT id FROM {0} " +
                     "WHERE id NOT IN ( " +
                     "SELECT id " +
                     "FROM genome_lists)").format(temp_table_name)

            self.cur.execute(query)

            missing_list_ids = []
            for (list_id,) in self.cur:
                missing_list_ids.append(list_id)

            if missing_list_ids:
                raise GenomeDatabaseError(
                    "Unknown genome list id(s) given. %s" % str(missing_list_ids))

            # Find any genome list ids that we dont have permission to view

            self.cur.execute("SELECT genome_id " +
                             "FROM genome_list_contents " +
                             "WHERE list_id in %s", (tuple(genome_list_ids),))

        except GenomeDatabaseError as e:
            raise e

        return [genome_id for (genome_id,) in self.cur.fetchall()]

    def deleteGenomeList(self, genome_list_ids):
        """Delete genome list and associated genomes.

        Parameters
        ----------
        genome_list_ids : iterable
            Unique identifier of genome lists in database.

        Returns
        -------
        bool
            True if successful.
        """

        for genome_list_id in genome_list_ids:
            try:
                edit_permission = self.permissionToModify(genome_list_id)
                if edit_permission is False:
                    raise GenomeDatabaseError(
                        "Insufficient permissions to delete genome list. Offending list id: {0}".format(genome_list_id))

                if not self._confirm("Are you sure you want to delete {0} lists (this action cannot be undone)".format(len(genome_list_ids))):
                    raise GenomeDatabaseError("User aborted database action.")

                list_genomes_ids = self.getGenomeIdsFromGenomeListIds(
                    [genome_list_id])
                self.editGenomeList(genome_list_id, list_genomes_ids, 'remove')
            except GenomeDatabaseError as e:
                raise e

        return True

    def printGenomeListsDetails(self, genome_list_ids):
        """Print genome list details.

        Parameters
        ----------
        genome_list_ids : iterable
            Unique identifier of genome lists in database.

        Returns
        -------
        list
            Column headers.
        list
            Content for each row.
        """

        try:
            if not genome_list_ids:
                raise GenomeDatabaseError(
                    "Unable to print genome details: No genomes given.")

            self.cur.execute(
                "SELECT lists.id, lists.name, lists.owned_by_root, users.username, count(contents.list_id) " +
                "FROM genome_lists as lists " +
                "LEFT OUTER JOIN users ON lists.owner_id = users.id " +
                "JOIN genome_list_contents as contents ON contents.list_id = lists.id " +
                "WHERE lists.id in %s " +
                "GROUP by lists.id, users.username " +
                "ORDER by lists.display_order asc, lists.id", (tuple(
                    genome_list_ids),)
            )

            header = ("list_id", "name", "owner", "genome_count")

            rows = []
            for (list_id, name, owned_by_root, username, genome_count) in self.cur:
                rows.append(
                    (list_id, name, ("root" if owned_by_root else username), genome_count))

        except GenomeDatabaseError as e:
            raise e

        return header, rows

    def permissionToModify(self, genome_list_id):
        """Check if user has permission to modify genome list.

        Parameters
        ----------
        genome_list_id : int
            Unique identifier of genome list in database.

        Returns
        -------
        bool
            True if has permission, else False.
        """

        try:
            self.cur.execute("SELECT owner_id, owned_by_root " +
                             "FROM genome_lists " +
                             "WHERE id = %s ", (genome_list_id,))

            result = self.cur.fetchone()

            if not result:
                raise GenomeDatabaseError(
                    "No genome list with id: %s" % str(genome_list_id))

            (owner_id, owned_by_root) = result

            if not self.currentUser.isRootUser():
                if owned_by_root or owner_id != self.currentUser.getUserId():
                    return False
            else:
                if not owned_by_root:
                    raise GenomeDatabaseError(
                        "Root user editing of other users lists not yet implemented.")

        except GenomeDatabaseError as e:
            raise e

        return True
