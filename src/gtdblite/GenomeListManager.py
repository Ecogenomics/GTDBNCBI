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

    # Creates a new genome list in the database
    #
    # Parameters:
    #     cur -
    #     genome_id_list - A list of genome ids to add to the new list.
    #     name - The name of the newly created list.
    #     description - A description of the newly created list.
    #     owner_id - The id of the user who will own this list.
    #     private - Bool that denotes whether this list is public or private (or none for not assigned).
    #
    # Returns:
    #    The genome list id of the newly created list.
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

    # True on success, false on failure/error.
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
                    "Insufficient pesrmissions to edit this genome list. Offending list id: %s" % genome_list_id)

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
                self.cur.execute("UPDATE genome_lists SET " + update_query +
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
                else:
                    raise GenomeDatabaseError(
                        "Unknown genome list edit operation: %s" % operation)
        except GenomeDatabaseError as e:
            raise e

        return True

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
