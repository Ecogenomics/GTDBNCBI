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

import hashlib
import logging
from multiprocessing import Pool

import prettytable

import Config
import Tools
from User import User
from GenomeDatabaseConnection import GenomeDatabaseConnection
from GenomeManager import GenomeManager
from GenomeListManager import GenomeListManager
from MetadataManager import MetadataManager
from GenomeFilter import GenomeFilter
from Exceptions import GenomeDatabaseError
from AlignedMarkerManager import AlignedMarkerManager


class GenomeDatabase(object):

    def __init__(self, threads=1, tab_table=False):
        self.logger = logging.getLogger()

        self.conn = GenomeDatabaseConnection()
        self.currentUser = None
        self.errorMessages = []
        self.warningMessages = []
        self.debugMode = False

        self.threads = threads
        self.pool = Pool(threads)

        self.tab_table = tab_table

        self.genomeCopyUserDir = None
        if Config.GTDB_GENOME_USR_DIR:
            self.genomeCopyUserDir = Config.GTDB_GENOME_USR_DIR

        self.defaultMarkerDatabaseName = 'user'

    #
    # Group: General Functions
    #
    # Function: ReportError
    # Sets the last error message of the database.
    #
    # Parameters:
    #     msg - The message to set.
    #
    # Returns:
    #   No return value.
    def ReportError(self, msg):
        self.errorMessages.append(str(msg))

    def GetErrors(self):
        return self.errorMessages

    def ClearErrors(self):
        self.errorMessages = []

    def ReportWarning(self, msg):
        self.warningMessages.append(str(msg))

    def GetWarnings(self):
        return self.warningMessages

    def ClearWarnings(self):
        self.warningMessages = []

    def PrintTable(self, header, rows):
        """Print table.

        Prints table information as either
        a tab-separated columns or as a
        prettytable.

        Parameters
        ----------
        header : iterable
            Column headers.
        rows : iterable
            Values for each column header.
        """

        if self.tab_table:
            print '\t'.join(header)
            for r in rows:
                print '\t'.join(map(str, r))
        else:
            table = prettytable.PrettyTable(header)
            table.align = 'l'
            table.hrules = prettytable.FRAME
            table.vrules = prettytable.NONE

            for r in rows:
                table.add_row(r)
            print table.get_string()

    # Function: SetDebugMode
    # Sets the debug mode of the database (at the moment its either on (non-zero) or off (zero))
    #
    # Parameters:
    #     debug_mode - The debug mode to set.
    #
    # Returns:
    #   No return value.
    def SetDebugMode(self, debug_mode):
        self.debugMode = debug_mode

    # TODO: This should not be here, techincally the backend is agnostic so
    # shouldn't assume command line.
    def Confirm(self, msg):
        raw = raw_input(msg + " (y/N): ")
        if raw.upper() == "Y":
            return True
        return False

    # Function: UserLogin
    # Log a user into the database (make the user the current user of the database).
    #
    # Parameters:
    #     username - The username of the user to login
    #
    # Returns:
    # Returns a User calls object on success (and sets the GenomeDatabase
    # current user), False otherwise.
    def UserLogin(self, username):
        try:
            if not self.conn.IsPostgresConnectionActive():
                raise GenomeDatabaseError(
                    "Unable to establish database connection")

            cur = self.conn.cursor()

            cur.execute("SELECT users.id, user_roles.id, user_roles.name "
                        "FROM users, user_roles " +
                        "WHERE users.role_id = user_roles.id " +
                        "AND users.username = %s", (username,))

            result = cur.fetchone()

            if not result:
                raise GenomeDatabaseError("User not found: %s" % username)

            (user_id, role_id, rolename) = result
            self.currentUser = User.createUser(
                user_id, username, rolename, role_id)

            return self.currentUser

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # Function: RootLogin
    # Log a user into the database as a root user (make the root user the current user of the database). Check
    # if the current user has permissions to do this.
    #
    # Parameters:
    #     username - The username of the user to login
    #
    # Returns:
    # Returns a User calls object on success (and sets the GenomeDatabase
    # current user), False otherwise.
    def RootLogin(self, username):
        try:
            if not self.conn.IsPostgresConnectionActive():
                raise GenomeDatabaseError(
                    "Unable to establish database connection")

            cur = self.conn.cursor()
            query = "SELECT id, has_root_login FROM users WHERE username = %s"
            cur.execute(query, [username])
            result = cur.fetchone()
            cur.close()
            if result:
                (_userid, has_root_login) = result
                if not has_root_login:
                    raise GenomeDatabaseError(
                        "You do not have sufficient permissions to logon as the root user.")

                self.currentUser = User.createRootUser(username)
                return self.currentUser

            raise GenomeDatabaseError("User %s not found." % username)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # Function: AddUser

    # Add a new user to the database.
    #
    # Parameters:
    #     username - The username of the user to login
    #     usertype - The role of the new user
    #
    # Returns:
    #   True on success, False otherwise.
    def AddUser(self, username, rolename=None, has_root=False):
        try:
            if rolename is None:
                rolename = 'user'

            if (not self.currentUser.isRootUser()):
                if has_root:
                    raise GenomeDatabaseError(
                        "Only the root user may grant root access to new users.")

                if rolename == 'admin':
                    raise GenomeDatabaseError(
                        "Only the root user may create admin accounts.")

                if not(self.currentUser.getRolename() == 'admin' and rolename == 'user'):
                    raise GenomeDatabaseError(
                        "Only admins (and root) can create user accounts.")

            cur = self.conn.cursor()

            cur.execute(
                "SELECT username from users where username = %s", (username,))

            if len(cur.fetchall()) > 0:
                raise GenomeDatabaseError(
                    "User %s already exists in the database." % username)

            cur.execute("INSERT into users (username, role_id, has_root_login) (" +
                        "SELECT %s, id, %s " +
                        "FROM user_roles " +
                        "WHERE name = %s)", (username, has_root, rolename))

            self.conn.commit()
            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return False
        except:
            self.conn.rollback()
            raise

    def EditUser(self, username, rolename=None, has_root=None):
        try:
            cur = self.conn.cursor()

            if (not self.currentUser.isRootUser()):
                raise GenomeDatabaseError(
                    "Only the root user may edit existing accounts.")

                # The following may be useful in the future if roles change, but at the moment,
                # only the root user can make any meaningful user edits
                """
                if has_root is not None:
                    raise GenomeDatabaseError("Only the root user may edit the root access of users.")

                if rolename == 'admin':
                    raise GenomeDatabaseError("Only the root user may create admin accounts.")

                cur.execute("SELECT users.id, user_roles.id, user_roles.name "
                    "FROM users, user_roles " +
                    "WHERE users.role_id = user_roles.id " +
                    "AND users.username = %s", (username, ))

                result = cur.fetchone()

                if not result:
                    raise GenomeDatabaseError("User not found: %s" % username)

                (user_id, current_role_id, current_rolename) = result
                if current_rolename == 'admin':
                    raise GenomeDatabaseError("Only the root user may edit current admin accounts.")
                """

            conditional_queries = []
            params = []

            if rolename is not None:
                conditional_queries.append(
                    " role_id = (SELECT id from user_roles where name = %s) ")
                params.append(rolename)

            if has_root is not None:
                conditional_queries.append(" has_root_login = %s ")
                params.append(has_root)

            if params:
                cur.execute("UPDATE users " +
                            "SET " + ','.join(conditional_queries) + " "
                            "WHERE username = %s", params + [username])

            self.conn.commit()
            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return False
        except:
            self.conn.rollback()
            raise

    # Function: GetUserIdFromUsername
    # Get a user id from a given username.
    #
    # Parameters:
    #     username - The username of the user to get an id for.
    #
    # Returns:
    #     The id of the user if successful, False on failure.
    def GetUserIdFromUsername(self, username):
        try:
            cur = self.conn.cursor()
            cur.execute(
                "SELECT id FROM users WHERE username = %s", (username,))
            result = cur.fetchone()

            if not result:
                raise GenomeDatabaseError("Username not found.")

            (user_id,) = result
            return user_id
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False
    #
    # Group: User Permission Functions
    #

    # Function: isCurrentUserRoleHigherThanUser
    # Checks if the current user is a higher user type than the specified user.
    #
    # Parameters:
    #     user_id - The id of the user to compare user types to.
    #
    # Returns:
    # True if the current user is a high user type than the user specified.
    # False otherwise. None on error.
    def isCurrentUserRoleHigherThanUser(self, user_id):
        try:
            cur = self.conn.cursor()
            cur.execute("SELECT type_id FROM users WHERE id = %s", (user_id,))
            result = cur.fetchone()

            if not result:
                raise GenomeDatabaseError("User not found.")

            (type_id,) = result
            if self.isRootUser() or (self.currentUser.getTypeId() < type_id):
                return True

            return False

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None

    def AddGenomes(self, batchfile, checkm_file, modify_genome_list_id=None,
                   new_genome_list_name=None):
        """Add genomes to database.

        Parameters
        ----------
        batchfile : str
            Name of file describing genomes to add.
        checkm_file : str
            Name of file containing CheckM quality estimates.

        Returns
        -------
        bool
            True if genomes added without error.
        """
        try:
            self.logger.info('Adding genomes to database.')

            cur = self.conn.cursor()

            if modify_genome_list_id is not None:
                if new_genome_list_name is not None:
                    raise GenomeDatabaseError(
                        "Unable to both modify and create genome lists at the same time.")

                genome_list_mngr = GenomeListManager(cur, self.currentUser)
                has_permission = genome_list_mngr.permissionToModify(
                    modify_genome_list_id)
                if has_permission is None:
                    raise GenomeDatabaseError(
                        "Unable to add genomes to list %s." % modify_genome_list_id)
                elif not has_permission:
                    raise GenomeDatabaseError(
                        "Insufficient permissions to add genomes to list %s." % modify_genome_list_id)

            if new_genome_list_name is not None:
                owner_id = None
                if not self.currentUser.isRootUser():
                    owner_id = self.currentUser.getUserId()

                genome_list_mngr = GenomeListManager(cur, self.currentUser)
                modify_genome_list_id = genome_list_mngr.addGenomeList([],
                                                                       new_genome_list_name,
                                                                       "",
                                                                       owner_id,
                                                                       True)
                if modify_genome_list_id is None:
                    raise GenomeDatabaseError(
                        "Unable to create the new genome list.")

            # add genomes to database
            genome_mngr = GenomeManager(cur, self.currentUser, self.threads)
            genome_ids = genome_mngr.addGenomes(checkm_file, batchfile)

            if modify_genome_list_id is not None:
                genome_list_mngr = GenomeListManager(cur, self.currentUser)
                bSuccess = genome_list_mngr.editGenomeList(modify_genome_list_id,
                                                           genome_ids=genome_ids,
                                                           operation='add')
                if not bSuccess:
                    raise GenomeDatabaseError(
                        "Unable to add genomes to genome list.")

            # all genomes were process successfully so move them into the GTDB
            # directory structure
            self.logger.info("Moving files to GTDB directory structure.")
            genome_mngr.moveGenomes(genome_ids)

            self.conn.commit()
            self.logger.info('Done.')
            return True
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False
        except:
            self.conn.rollback()
            raise

    def list_markers(self, cur=None, library=None):
        cur.execute("SELECT id_in_database FROM markers m " +
                    "LEFT OUTER JOIN marker_databases md ON m.marker_database_id = md.id " +
                    "WHERE md.name like '{0}'".format(library))
        listmarkers = [hmm_id for (hmm_id,) in cur.fetchall()]
        return listmarkers

    # True if has permission. False if doesn't. None on error.
    def DeleteGenomes(self, batchfile=None, external_ids=None, reason=None):
        try:
            if external_ids is None:
                external_ids = []

            if batchfile:
                fh = open(batchfile, "rb")
                external_ids.extend(
                    [line.rstrip().split('\t')[0] for line in fh])
                fh.close()

            genome_ids = self.ExternalGenomeIdsToGenomeIds(external_ids)

            cur = self.conn.cursor()
            genome_mngr = GenomeManager(cur, self.currentUser, 1)
            genome_mngr.deleteGenomes(batchfile, genome_ids, reason)

            self.conn.commit()
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

    def GetAllUserGenomeIds(self):
        """Get genome identifiers for all user genomes.

        Returns
        -------
        list
            Database identifiers for all user genomes.
        """

        try:
            cur = self.conn.cursor()

            cur.execute("SELECT id " +
                        "FROM genome_sources " +
                        "WHERE name = %s", (self.defaultMarkerDatabaseName,))

            user_source_id = cur.fetchone()[0]

            cur.execute("SELECT id " +
                        "FROM genomes " +
                        "WHERE genome_source_id = %s", (user_source_id,))

            result_ids = []
            for (genome_id,) in cur:
                result_ids.append(genome_id)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return result_ids

    def GetAllGenomeIds(self):
        """Get genome identifiers for all genomes.

        Returns
        -------
        list
            Database identifiers for all genomes.
        """

        try:
            cur = self.conn.cursor()

            query = "SELECT id FROM genomes"
            cur.execute(query)

            result_ids = []
            for (genome_id,) in cur:
                result_ids.append(genome_id)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return result_ids

    # True on success. False on failure/error.
    def ViewGenomes(self, batchfile=None, external_ids=None):
        try:
            genome_ids = []
            if external_ids is None and batchfile is None:
                genome_ids = self.GetAllGenomeIds()
            else:
                if external_ids is None:
                    external_ids = []
                if batchfile:
                    try:
                        fh = open(batchfile, "rb")
                    except:
                        raise GenomeDatabaseError(
                            "Cannot open batchfile: " + batchfile)

                    for line in fh:
                        line = line.rstrip()
                        external_ids.append(line)

                genome_ids = self.ExternalGenomeIdsToGenomeIds(external_ids)
                if genome_ids is None:
                    raise GenomeDatabaseError("Can not retrieve genome ids.")

            return self.PrintGenomesDetails(genome_ids)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # List of genome ids on success. False on error.
    def ExternalGenomeIdsToGenomeIds(self, external_ids):
        '''
        Function
        Check if All external ids are stored in GTDB

        :param external_ids: List of ids to delete
        '''
        try:
            cur = self.conn.cursor()

            map_sources_to_ids = {}

            for external_id in external_ids:
                try:
                    (source_prefix, id_at_source) = external_id.split("_", 1)
                except ValueError:
                    raise GenomeDatabaseError(
                        "All genome ids must have the form <prefix>_<id>. Offending id: %s" % str(external_id))

                if source_prefix not in map_sources_to_ids:
                    map_sources_to_ids[source_prefix] = {}
                map_sources_to_ids[source_prefix][id_at_source] = external_id

            temp_table_name = Tools.generateTempTableName()

            if len(map_sources_to_ids.keys()):
                cur.execute("CREATE TEMP TABLE %s (prefix text)" %
                            (temp_table_name,))
                query = "INSERT INTO {0} (prefix) VALUES (%s)".format(
                    temp_table_name)
                cur.executemany(query, [(x,)
                                        for x in map_sources_to_ids.keys()])
            else:
                raise GenomeDatabaseError(
                    "No genome sources found for these ids. %s" % str(external_ids))

            # Find any given tree prefixes that arent in the genome sources
            query = ("SELECT prefix FROM {0} " +
                     "WHERE prefix NOT IN ( " +
                     "SELECT external_id_prefix " +
                     "FROM genome_sources)").format(temp_table_name)

            cur.execute(query)

            missing_genome_sources = {}
            for (query_prefix,) in cur:
                missing_genome_sources[query_prefix] = map_sources_to_ids[
                    query_prefix].values()

            if len(missing_genome_sources.keys()):
                errors = []
                for (source_prefix, offending_ids) in missing_genome_sources.items():
                    errors.append("(%s) %s" %
                                  (source_prefix, str(offending_ids)))
                raise GenomeDatabaseError("Cannot find the relevant genome source id for the following ids, check the IDs are correct: " +
                                          ", ".join(errors))

            # All genome sources should be good, find ids
            result_ids = []
            for source_prefix in map_sources_to_ids.keys():

                # Create a table of requested external ids from this genome
                # source
                temp_table_name = Tools.generateTempTableName()
                cur.execute(
                    "CREATE TEMP TABLE %s (id_at_source text)" % (temp_table_name,))
                query = "INSERT INTO {0} (id_at_source) VALUES (%s)".format(
                    temp_table_name)
                cur.executemany(
                    query, [(x,) for x in map_sources_to_ids[source_prefix].keys()])

                # Check to see if there are any that don't exist
                query = ("SELECT id_at_source FROM {0} " +
                         "WHERE id_at_source NOT IN ( " +
                         "SELECT id_at_source " +
                         "FROM genomes, genome_sources " +
                         "WHERE genome_source_id = genome_sources.id " +
                         "AND external_id_prefix = %s)").format(temp_table_name)

                cur.execute(query, (source_prefix,))

                missing_ids = []
                for (id_at_source,) in cur:
                    missing_ids.append(source_prefix + "_" + id_at_source)

                if missing_ids:
                    raise GenomeDatabaseError(
                        "Cannot find the the following genome ids, check the IDs are correct: %s" % str(missing_ids))

                # All exist, so get their ids.
                query = ("SELECT genomes.id FROM genomes, genome_sources " +
                         "WHERE genome_source_id = genome_sources.id " +
                         "AND id_at_source IN ( " +
                         "SELECT id_at_source " +
                         "FROM {0} )" +
                         "AND external_id_prefix = %s").format(temp_table_name)

                cur.execute(query, (source_prefix,))

                for (genome_id,) in cur:
                    result_ids.append(genome_id)

            return result_ids

        except GenomeDatabaseError as e:
            raise e

    def PrintGenomesDetails(self, genome_id_list):
        """Print genome details.

        Parameters
        ----------
        genome_id_list : iterable
            Unique identifier of genomes in database.

        Returns
        -------
        bool
            True if successful, else False.
        """

        try:
            if not genome_id_list:
                raise GenomeDatabaseError(
                    "Unable to print genomes. No genomes found.")

            cur = self.conn.cursor()

            columns = "genomes.id, genomes.name, description, owned_by_root, username, " + \
                "external_id_prefix || '_' || id_at_source as external_id, date_added"

            cur.execute("SELECT " + columns + " FROM genomes " +
                        "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
                        "JOIN genome_sources AS sources ON genome_source_id = sources.id " +
                        "AND genomes.id in %s " +
                        "ORDER BY genomes.id ASC", (tuple(genome_id_list),))

            # print table
            header = (
                "genome_id", "name", "description", "owner", "data_added")

            rows = []
            for (_genome_id, name, description, owned_by_root, username, external_id, date_added) in cur:
                rows.append((external_id, name, description,
                             ("root" if owned_by_root else username),
                             date_added.date()))

            self.PrintTable(header, rows)
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

    def AddMarkerWorking(self, cur, marker_file_path, name, desc, marker_set_id=None, force_overwrite=False,
                         database=None, id_in_database=None):
        try:
            try:
                marker_fh = open(marker_file_path, "rb")
            except:
                raise GenomeDatabaseError(
                    "Cannot open Marker file: " + marker_file_path)

            seen_name_line = False
            model_length = None

            m = hashlib.sha256()
            for line in marker_fh:
                if line[:4] == 'NAME':
                    if seen_name_line:
                        raise GenomeDatabaseError(
                            "Marker file contains more than one model. Offending file: " + marker_file_path)
                    seen_name_line = True
                elif line[:4] == 'LENG':
                    try:
                        model_length = int(line[4:])
                    except:
                        raise GenomeDatabaseError(
                            "Unable to convert model length into integer value. Offending line: %s. Offending file %s." % (line, marker_file_path))
                m.update(line)

            if model_length is None:
                raise GenomeDatabaseError(
                    "Model file does not give specify marker length. Offending file %s." % marker_file_path)

            if model_length <= 0:
                raise GenomeDatabaseError("Model file specifies invalid marker length. Length: %i. Offending file %s." % (
                    model_length, marker_file_path))

            marker_sha256_checksum = m.hexdigest()
            marker_fh.close()

            if database is None:
                database = self.defaultMarkerDatabaseName

            if marker_set_id is not None:
                if self.GetMarkerIdListFromMarkerSetId(marker_set_id) is False:
                    raise GenomeDatabaseError(
                        "Unable to add marker to set %s." % marker_set_id)

            if marker_set_id is not None:
                has_permission = self.HasPermissionToEditMarkerSet(
                    marker_set_id)
                if has_permission is None:
                    raise GenomeDatabaseError(
                        "Unable to add marker to set %s." % marker_set_id)
                elif not has_permission:
                    raise GenomeDatabaseError(
                        "Insufficient permission to add marker to marker set %s." % marker_set_id)

            cur.execute(
                "SELECT id, external_id_prefix, user_editable FROM marker_databases WHERE name = %s", (database,))
            database_id = None

            for (this_database_id, _external_id_prefix, user_editable) in cur:
                if (not user_editable):
                    if id_in_database is None:
                        raise GenomeDatabaseError(
                            "Cannot auto generate ids in databases for the %s marker database." % database)
                    if (not self.currentUser.isRootUser()):
                        raise GenomeDatabaseError(
                            "Only the root user can add markers to the %s marker database." % database)
                database_id = this_database_id
                break

            if database_id is None:
                raise GenomeDatabaseError(
                    "Could not find the %s marker database." % database)

            if id_in_database is None:
                cur.execute(
                    "SELECT id_in_database FROM markers WHERE marker_database_id = %s order by id_in_database::int desc", (database_id,))
                last_id = None
                for (last_id_in_database,) in cur:
                    last_id = last_id_in_database
                    break

                cur.execute(
                    "SELECT last_auto_id FROM marker_databases WHERE id = %s ", (database_id,))
                for (last_auto_id,) in cur:
                    if last_id is None:
                        last_id = last_auto_id
                    else:
                        last_id = max(int(last_id), int(last_auto_id))
                    break

                # Generate a new id (for user editable lists only)
                if (last_id is None):
                    new_id = 1
                else:
                    new_id = int(last_id) + 1

                if id_in_database is None:
                    id_in_database = str(new_id)

                cur.execute(
                    "UPDATE marker_databases set last_auto_id = %s where id = %s", (new_id, database_id))

            owner_id = None
            if not self.currentUser.isRootUser():
                owner_id = self.currentUser.getUserId()

            cur.execute(
                "SELECT id FROM markers WHERE marker_database_id = %s AND id_in_database = %s", (database_id, id_in_database))

            result = cur.fetchall()

            columns = "(name, description, owned_by_root, owner_id, marker_file_location, " + \
                "marker_file_sha256, marker_database_id, id_in_database, size)"

            if len(result):
                if force_overwrite:
                    raise GenomeDatabaseError(
                        "Force overwrite not implemented yet")
                else:
                    raise GenomeDatabaseError(
                        "Marker database '%s' already contains id '%s'. Use -f to force an overwrite." % (database, id_in_database))

            cur.execute("INSERT INTO markers " + columns + " "
                        "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s) " +
                        "RETURNING id",
                        (name, desc, self.currentUser.isRootUser(), owner_id, marker_file_path, marker_sha256_checksum, database_id, id_in_database, model_length))

            (marker_id,) = cur.fetchone()

            # TODO: Add to marker set if needed

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return marker_id

    def GetAllMarkerIds(self):
        try:
            cur = self.conn.cursor()

            query = "SELECT id FROM markers"
            cur.execute(query)

            result_ids = []
            for (marker_id,) in cur:
                result_ids.append(marker_id)

            return result_ids

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def ExternalMarkerIdsToMarkerIds(self, external_ids):
        try:
            cur = self.conn.cursor()

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
                cur.execute("CREATE TEMP TABLE %s (prefix text)" %
                            (temp_table_name,))
                query = "INSERT INTO {0} (prefix) VALUES (%s)".format(
                    temp_table_name)
                cur.executemany(query, [(x,)
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

            cur.execute(query)

            missing_marker_sources = {}
            for (query_prefix,) in cur:
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
                cur.execute(
                    "CREATE TEMP TABLE %s (id_in_database text)" % (temp_table_name,))
                query = "INSERT INTO {0} (id_in_database) VALUES (%s)".format(
                    temp_table_name)
                cur.executemany(
                    query, [(x,) for x in map_databases_to_ids[database_prefix].keys()])

                # Check to see if there are any that don't exist
                query = ("SELECT id_in_database FROM {0} " +
                         "WHERE id_in_database NOT IN ( " +
                         "SELECT id_in_database " +
                         "FROM markers, marker_databases " +
                         "WHERE marker_database_id = marker_databases.id " +
                         "AND external_id_prefix = %s)").format(temp_table_name)

                cur.execute(query, (database_prefix,))

                missing_ids = []
                for (id_in_database,) in cur:
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

                cur.execute(query, (database_prefix,))

                for (marker_id,) in cur:
                    result_ids.append(marker_id)

            return result_ids

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def ViewMarkers(self, batchfile=None, external_ids=None):
        try:
            marker_ids = []
            if external_ids is None and batchfile is None:
                marker_ids = self.GetAllMarkerIds()
            else:
                if external_ids is None:
                    external_ids = []
                if batchfile:
                    try:
                        fh = open(batchfile, "rb")
                    except:
                        raise GenomeDatabaseError(
                            "Cannot open batchfile: " + batchfile)

                    for line in fh:
                        line = line.rstrip()
                        external_ids.append(line)

                marker_ids = self.ExternalMarkerIdsToMarkerIds(external_ids)
                if marker_ids is False:
                    raise GenomeDatabaseError("Can not retrieve marker ids.")

            return self.PrintMarkerDetails(marker_ids)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def PrintMarkerDetails(self, marker_id_list):
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

            cur = self.conn.cursor()

            columns = "markers.id, markers.name, description, " + \
                "external_id_prefix || '_' || id_in_database as external_id, size"

            cur.execute("SELECT " + columns + " FROM markers " +
                        "LEFT OUTER JOIN users ON markers.owner_id = users.id " +
                        "JOIN marker_databases AS databases ON marker_database_id = databases.id " +
                        "AND markers.id in %s " +
                        "ORDER BY markers.id ASC", (tuple(marker_id_list),))

            # print table
            header = ("marker_id", "name", "description", "size (nt)")

            rows = []
            for (_marker_id, name, description, external_id, size) in cur:
                rows.append((external_id, name, description, size))

            self.PrintTable(header, rows)

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # Function: GetMarkerIdListFromMarkerSetId
    # Given a marker set id, return all the ids of the markers contained within that marker set.
    #
    # Parameters:
    #     marker_set_id - The marker set id of the marker set whose contents needs to be retrieved.
    #
    # Returns:
    # A list of all the marker ids contained within the specified marker set,
    # False on failure.
    def GetMarkerIdListFromMarkerSetId(self, marker_set_id):

        cur = self.conn.cursor()

        cur.execute("SELECT id, owner_id, owned_by_root, private " +
                    "FROM marker_sets " +
                    "WHERE id = %s ", (tuple(marker_set_id),))

        result = cur.fetchone()

        if not result:
            self.ReportError("No marker set with id: %s" % str(marker_set_id))
            return None
        else:
            (_list_id, owner_id, owned_by_root, private) = result
            if private and (not self.currentUser.isRootUser()) and (owned_by_root or owner_id != self.currentUser.getUserId()):
                self.ReportError(
                    "Insufficient permission to view marker set: %s" % str(marker_set_id))
                return None

        cur.execute("SELECT marker_id " +
                    "FROM marker_set_contents " +
                    "WHERE set_id = %s ", (tuple(marker_set_id),))

        return [marker_id for (marker_id,) in cur.fetchall()]

    def FindUncalculatedMarkersForGenomeId(self, genome_id, marker_ids):

        cur = self.conn.cursor()

        cur.execute("SELECT marker_id, sequence " +
                    "FROM aligned_markers " +
                    "WHERE genome_id = %s ", (genome_id,))

        marker_id_dict = dict(cur.fetchall())

        return [x for x in marker_ids if x not in marker_id_dict]

    def MakeTreeData(self, marker_ids, genome_ids,
                     directory, prefix,
                     comp_threshold, cont_threshold,
                     taxa_filter,
                     guaranteed_genome_list_ids, guaranteed_genome_ids,
                     alignment, individual,
                     build_tree=True):

        self.logger.info('Making tree for %d genomes using %d marker genes.' % (
            len(genome_ids), len(marker_ids)))

        try:
            gf = GenomeFilter()
            genomes_to_retain, chosen_markers_order, chosen_markers = gf.FilterTreeData(self, marker_ids, genome_ids,
                                                                                        comp_threshold, cont_threshold,
                                                                                        taxa_filter,
                                                                                        guaranteed_genome_list_ids, guaranteed_genome_ids,
                                                                                        directory)

            aligned_mngr = AlignedMarkerManager(self.threads)
            aligned_mngr.calculateAlignedMarkerSets(
                genomes_to_retain, marker_ids)

            gf.writeTreeFiles(
                self, marker_ids, genomes_to_retain, directory, prefix, chosen_markers_order, chosen_markers, alignment, individual)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        if build_tree:
            pass
            # To Do: should infer tree here

        self.logger.info('Done.')

        return True

    def GetAlignedMarkersCountForGenomes(self, genome_ids, marker_ids):

        cur = self.conn.cursor()

        cur.execute("SELECT genome_id, count(marker_id) " +
                    "FROM aligned_markers " +
                    "WHERE genome_id in %s " +
                    "AND marker_id in %s " +
                    "AND sequence IS NOT NULL " +
                    "GROUP BY genome_id", (tuple(genome_ids), tuple(marker_ids)))

        return dict(cur.fetchall())

    def CreateGenomeList(self, batchfile, external_ids, name, description, private=None):
        try:
            cur = self.conn.cursor()

            genome_id_list = []

            if not external_ids:
                external_ids = []

            if batchfile:
                fh = open(batchfile, "rb")
                for line in fh:
                    if line[0] != '#':
                        external_id = line.rstrip().split('\t')[0]
                        external_ids.append(external_id)

            if not external_ids:
                raise GenomeDatabaseError(
                    "No genomes provided to create a genome list.")

            genome_id_list = self.ExternalGenomeIdsToGenomeIds(external_ids)
            if genome_id_list is False:
                raise GenomeDatabaseError(
                    "Unable to retreive genome ids for provided genomes.")

            owner_id = None
            if not self.currentUser.isRootUser():
                owner_id = self.currentUser.getUserId()

            genome_list_mngr = GenomeListManager(cur, self.currentUser)
            genome_list_id = genome_list_mngr.addGenomeList(genome_id_list,
                                                            name,
                                                            description,
                                                            owner_id,
                                                            private)
            if genome_list_id is False:
                raise GenomeDatabaseError("Unable to create new genome list.")

            self.conn.commit()

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return genome_list_id

    # Function: GetVisibleGenomeLists
    # Get all the genome lists that the current user can see.
    #
    # Parameters:
    #     owner_id - Get visible genome lists owned by this user with this id. If not specified, get all root owned lists.
    #     include_private - If true, private genome are also shown.
    #
    # Returns:
    #   A list containing a tuple for each visible genome list. The tuple contains the genome list id, genome list name, genome list description,
    # and username of the owner of the list (id, name, description, username).

    def GetVisibleGenomeListsByOwner(self, owner_id=None, include_private=False):
        """
        Get all genome list owned by owner_id which the current user is allowed
        to see. If owner_id is None, return all visible genome lists for the
        current user.
        """
        cur = self.conn.cursor()

        conditional_query = ""
        params = []

        if owner_id is None:
            conditional_query += "AND owned_by_root is True "
        else:
            conditional_query += "AND owner_id = %s "
            params.append(owner_id)

        if not self.currentUser.isRootUser():
            privacy_condition = "private = False"
            if include_private:
                privacy_condition = "(private = False OR private = True)"

            conditional_query += "AND (" + \
                privacy_condition + " OR owner_id = %s)"
            params.append(self.currentUser.getUserId())

        cur.execute("SELECT id " +
                    "FROM genome_lists " +
                    "WHERE 1 = 1 " +
                    conditional_query, params)

        return [list_id for (list_id,) in cur]

    def GetAllVisibleGenomeListIds(self, include_private=False):
        cur = self.conn.cursor()

        conditional_query = ""
        params = []
        if not self.currentUser.isRootUser():
            privacy_condition = "private = False"
            if include_private:
                privacy_condition = "(private = False OR private = True)"

            conditional_query += "AND (" + \
                privacy_condition + " OR owner_id = %s)"
            params.append(self.currentUser.getUserId())

        cur.execute("SELECT id " +
                    "FROM genome_lists " +
                    "WHERE 1 = 1 " +
                    conditional_query, params)

        return [list_id for (list_id,) in cur]

    def ViewGenomeListsContents(self, list_ids):
        try:
            cur = self.conn.cursor()
            genome_list_mngr = GenomeListManager(cur, self.currentUser)
            header, rows = genome_list_mngr.viewGenomeListsContents(list_ids)
            if not header:
                raise GenomeDatabaseError(
                    "Unable to view genome lists")
            else:
                self.PrintTable(header, rows)

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # True on success, false on failure.
    def EditGenomeList(self, genome_list_id, batchfile=None, genomes_external_ids=None, operation=None, name=None, description=None, private=None):
        """Edit an existing genome list in the database.

        Parameters
        ----------
        genome_list_id : int
            Identifier of genome list in database.
        batchfile : str
            Filename of batch file describing genomes to modify.
        genomes_external_ids : list
            List of genomes to modify.
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
            cur = self.conn.cursor()

            if batchfile:
                if genomes_external_ids is None:
                    genomes_external_ids = []
                fh = open(batchfile, 'rb')
                for line in fh:
                    line = line.rstrip()
                    genomes_external_ids.append(line)
                fh.close()

            genome_ids = []
            if genomes_external_ids is not None:
                genome_ids = self.ExternalGenomeIdsToGenomeIds(
                    genomes_external_ids)

            if genome_ids is False:
                raise GenomeDatabaseError(
                    "Unable to retrive information for all genome ids.")

            genome_list_mngr = GenomeListManager(cur, self.currentUser)
            if not genome_list_mngr.editGenomeList(genome_list_id, genome_ids, operation, name, description, private):
                raise GenomeDatabaseError(
                    "Unable to edit genome list: %s" % genome_list_id)

            self.conn.commit()
        except GenomeDatabaseError as e:
            self.conn.rollback()
            self.ReportError(e.message)
            return False

        return True

    def deleteGenomeLists(self, list_ids=None):
        '''
        Delete a list of genome lists

        :param list_ids: List of list IDs to delete 
        '''
        try:
            cur = self.conn.cursor()
            genome_list_mngr = GenomeListManager(cur, self.currentUser)
            if not genome_list_mngr.deleteGenomeList(list_ids):
                raise GenomeDatabaseError(
                    "Unable to edit genome lists")

            self.conn.commit()
        except GenomeDatabaseError as e:
            self.conn.rollback()
            self.ReportError(e.message)
            return False
        return True

    def CreateMarkerSet(self, batchfile, external_ids, name, description, private=None):
        try:
            cur = self.conn.cursor()

            marker_id_list = []

            if not external_ids:
                external_ids = []

            if batchfile:
                fh = open(batchfile, "rb")
                for line in fh:
                    if line[0] == '#':
                        continue
                    marker_id = line.rstrip().split('\t')[0]
                    external_ids.append(marker_id)

            if not external_ids:
                raise GenomeDatabaseError(
                    "No markers provided to create a marker set.")

            marker_id_list = self.ExternalMarkerIdsToMarkerIds(external_ids)
            if marker_id_list is False:
                raise GenomeDatabaseError(
                    "Unable to retreive marker ids for provided markers.")

            owner_id = None
            if not self.currentUser.isRootUser():
                owner_id = self.currentUser.getUserId()

            marker_list_id = self.CreateMarkerSetWorking(
                cur, marker_id_list, name, description, owner_id, private)
            if marker_list_id is False:
                raise GenomeDatabaseError("Unable to create new marker set.")

            self.conn.commit()

            return marker_list_id

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # True on success, false on failure/error.
    def CreateMarkerSetWorking(self, cur, marker_id_list, name, description, owner_id=None, private=None):
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
            cur.execute(
                query, (name, description, owner_id is None, owner_id, private))
            (marker_set_id,) = cur.fetchone()

            query = "INSERT INTO marker_set_contents (set_id, marker_id) VALUES (%s, %s)"
            cur.executemany(query, [(marker_set_id, x)
                                    for x in marker_id_list])

            return marker_set_id

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # True on success, false on failure/error.
    def EditMarkerSet(self, marker_set_id, batchfile=None, marker_external_ids=None, operation=None, name=None, description=None, private=None):

        cur = self.conn.cursor()

        if batchfile:
            if marker_external_ids is None:
                marker_external_ids = []
            for line in batchfile:
                line = line.rstrip()
                marker_external_ids.append(line)

        if marker_external_ids is not None:
            marker_external_ids = self.ExternalMarkerIdsToMarkerIds(
                marker_external_ids)

        if not self.EditMarkerSetWorking(cur, marker_set_id, marker_external_ids, operation, name, description, private):
            self.conn.rollback()
            return False

        self.conn.commit()
        return True

    # True on success, false on failure/error.
    def EditMarkerSetWorking(self, cur, marker_set_id, marker_ids=None, operation=None, name=None, description=None, private=None):
        try:
            edit_permission = self.HasPermissionToEditMarkerSet(marker_set_id)
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
                cur.execute("UPDATE marker_sets SET " + update_query +
                            " WHERE id = %s", params + [marker_set_id])

            temp_table_name = Tools.generateTempTableName()

            if operation is not None:

                if len(marker_ids) == 0:
                    raise GenomeDatabaseError(
                        "No marker ids given to perform '%s' operation." % operation)

                cur.execute("CREATE TEMP TABLE %s (id integer)" %
                            (temp_table_name,))
                query = "INSERT INTO {0} (id) VALUES (%s)".format(
                    temp_table_name)
                cur.executemany(query, [(x,) for x in marker_ids])

                if operation == 'add':
                    query = ("INSERT INTO marker_set_contents (set_id, marker_id) " +
                             "SELECT %s, id FROM {0} " +
                             "WHERE id NOT IN ( " +
                             "SELECT marker_id " +
                             "FROM marker_set_contents " +
                             "WHERE set_id = %s)").format(temp_table_name)
                    cur.execute(query, (marker_set_id, marker_set_id))
                elif operation == 'remove':
                    query = ("DELETE FROM marker_set_contents " +
                             "WHERE set_id = %s " +
                             "AND marker_id IN ( " +
                             "SELECT id " +
                             "FROM {0})").format(temp_table_name)
                    cur.execute(query, [marker_set_id])
                else:
                    raise GenomeDatabaseError(
                        "Unknown marker set edit operation: %s" % operation)

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # True if has permission, False if not. None on error.
    def HasPermissionToEditMarkerSet(self, marker_set_id):
        try:
            cur = self.conn.cursor()

            cur.execute("SELECT owner_id, owned_by_root " +
                        "FROM marker_sets " +
                        "WHERE id = %s ", (marker_set_id,))

            result = cur.fetchone()

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

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return None

    # True on success, false on failure/error.
    def PrintMarkerSetsDetails(self, marker_set_ids):
        try:
            cur = self.conn.cursor()

            if not marker_set_ids:
                raise GenomeDatabaseError(
                    "Unable to print marker set details: No marker sets given.")

            if not self.currentUser.isRootUser():
                cur.execute("SELECT id " +
                            "FROM marker_sets as sets " +
                            "WHERE sets.private = True " +
                            "AND sets.id in %s " +
                            "AND (owned_by_root = True OR owner_id != %s)", (tuple(marker_set_ids), self.currentUser.getUserId()))

                unviewable_set_ids = [set_id for (set_id,) in cur]
                if unviewable_set_ids:
                    raise GenomeDatabaseError(
                        "Insufficient privileges to view marker sets: %s." % str(unviewable_set_ids))

            cur.execute(
                "SELECT sets.id, sets.name, sets.description, sets.owned_by_root, users.username, count(contents.set_id) " +
                "FROM marker_sets as sets " +
                "LEFT OUTER JOIN users ON sets.owner_id = users.id " +
                "JOIN marker_set_contents as contents ON contents.set_id = sets.id " +
                "WHERE sets.id in %s " +
                "GROUP by sets.id, users.username " +
                "ORDER by sets.id asc ", (tuple(marker_set_ids),)
            )

            # print table
            header = ("set_id", "name", "description", "owner", "marker_count")

            rows = []
            for (set_id, name, description, owned_by_root, username, marker_count) in cur:
                rows.append(
                    (set_id, name, description, ("root" if owned_by_root else username), marker_count))

            self.PrintTable(header, rows)

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def GetVisibleMarkerSetsByOwner(self, owner_id=None, include_private=False):
        """
        Get all marker sets owned by owner_id which the current user is allowed
        to see. If owner_id is None, return all visible marker sets for the
        current user.
        """

        cur = self.conn.cursor()

        conditional_query = ""
        params = []

        if owner_id is None:
            conditional_query += "AND owned_by_root is True "
        else:
            conditional_query += "AND owner_id = %s "
            params.append(owner_id)

        if not self.currentUser.isRootUser():
            privacy_condition = "private = False"
            if include_private:
                privacy_condition = "(private = False OR private = True)"

            conditional_query += "AND (" + \
                privacy_condition + " OR owner_id = %s)"
            params.append(self.currentUser.getUserId())

        cur.execute("SELECT id " +
                    "FROM marker_sets " +
                    "WHERE 1 = 1 " +
                    conditional_query, params)

        return [set_id for (set_id,) in cur]

    # Returns list of marker set id. False on failure/error.
    def GetAllVisibleMarkerSetIds(self, include_private=False):
        cur = self.conn.cursor()

        conditional_query = ""
        params = []

        if not self.currentUser.isRootUser():
            privacy_condition = "private = False"
            if include_private:
                privacy_condition = "(private = False OR private = True)"

            conditional_query += "AND (" + \
                privacy_condition + " OR owner_id = %s)"
            params.append(self.currentUser.getUserId())

        cur.execute("SELECT id " +
                    "FROM marker_sets " +
                    "WHERE 1 = 1 " +
                    conditional_query, params)

        return [set_id for (set_id,) in cur]

    # Returns list of marker set id. False on failure/error.
    def ViewMarkerSetsContents(self, marker_set_ids):
        try:
            marker_ids = self.GetMarkerIdListFromMarkerSetId(marker_set_ids)

            if marker_ids is None:
                raise GenomeDatabaseError(
                    "Unable to view marker set. Can not retrieve marker IDs for sets: %s" % str(marker_set_ids))

            if not self.PrintMarkerDetails(marker_ids):
                raise GenomeDatabaseError(
                    "Unable to view marker set. Printing to screen failed of marker ids. %s" % str(marker_ids))

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def viewMetadata(self):
        try:
            metaman = MetadataManager()
            metaman.viewMetadata()
            return True
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def exportMetadata(self, path):
        try:
            metaman = MetadataManager()
            metaman.exportMetadata(path)
            return True
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def importMetadata(self, table=None, field=None, typemeta=None, metafile=None):
        try:
            metaman = MetadataManager()
            metaman.importMetadata(table, field, typemeta, metafile)
            return True
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def createMetadata(self, path):
        try:
            metaman = MetadataManager()
            metaman.createMetadata(path)
            return True
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def reportStats(self):
        """Report general database statistics."""

        try:
            cur = self.conn.cursor()

            cur.execute("SELECT name, external_id_prefix, " +
                        "(SELECT COUNT(*) "
                        "FROM genomes "
                        "WHERE genomes.genome_source_id = genome_sources.id) "
                        "FROM genome_sources " +
                        "ORDER BY id")
            genome_counts = cur.fetchall()

            print '\t'.join(('Genome Source', 'Prefix', 'Genome Count'))
            for tup in genome_counts:
                print '\t'.join(map(str, list(tup)))

            print ''
            print 'Total genomes: %d' % (sum([x[2] for x in genome_counts]))

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True
