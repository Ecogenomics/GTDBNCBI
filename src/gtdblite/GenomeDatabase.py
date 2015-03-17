import hashlib
import shutil
import os
import datetime
import time
import random

from gtdblite.User import User
from gtdblite.GenomeDatabaseConnection import GenomeDatabaseConnection


class GenomeDatabaseError(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


class GenomeDatabase(object):
    def __init__(self):
        self.conn = GenomeDatabaseConnection()
        self.currentUser = None
        self.errorMessages = []
        self.warningMessages = []
        self.debugMode = False
        self.genomeCopyDir = None
        self.defaultGenomeSourceName = 'user'
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

    # Function: UserLogin
    # Log a user into the database (make the user the current user of the database).
    #
    # Parameters:
    #     username - The username of the user to login
    #     password - The password of the user
    #
    # Returns:
    #   Returns a User calls object on success (and sets the GenomeDatabase current user), None otherwise.
    def UserLogin(self, username):
        if not self.conn.IsPostgresConnectionActive():
            self.ReportError("Unable to establish database connection")
            return None

        cur = self.conn.cursor()
        query = "SELECT id, role_id FROM users WHERE username = %s"
        cur.execute(query, [username])
        result = cur.fetchone()
        cur.close()
        if result:
            (userid, role_id) = result
            self.currentUser = User.createUser(result[0], username, result[1])
            return self.currentUser
        else:
            self.ReportError("User not found: %s" % username)
        return None

    # Function: RootLogin
    # Log a user into the database as a root user (make the root user the current user of the database). Check
    # if the current user has permissions to do this.
    #
    # Parameters:
    #     username - The username of the user to login
    #     password - The password of the user
    #
    # Returns:
    #   Returns a User calls object on success (and sets the GenomeDatabase current user), None otherwise.
    def RootLogin(self, username):
        if not self.conn.IsPostgresConnectionActive():
            self.ReportError("Unable to establish database connection")
            return None

        cur = self.conn.cursor()
        query = "SELECT id, has_root_login FROM users WHERE username = %s"
        cur.execute(query, [username])
        result = cur.fetchone()
        cur.close()
        if result:
            (userid, has_root_login) = result
            if has_root_login:
                self.currentUser = User.createRootUser()
                return self.currentUser
            else:
                self.ReportError("You do not have sufficient permissions to logon as the root user.")
        else:
            self.ReportError("User %s not found." % username)
        return None

    # Function: CreateUser

    # Create a new user for the database.
    #
    # Parameters:
    #     username - The username of the user to login
    #     userTypeId - The id of the type of user to create
    #
    # Returns:
    #   True on success, False otherwise.
    def CreateUser(self, username, userTypeId):

        currentUser = self.currentUser

        if not self.conn.IsPostgresConnectionActive():
            self.ReportError("Unable to establish database connection")
            return False

        if not currentUser:
            self.ReportError("You need to be logged in to create a user")
            return False


        if (not currentUser.isRootUser()) or userTypeId <= self.currentUser.getTypeId():
            self.ReportError("Cannot create a user with same or higher level privileges")
            return False

        cur = self.conn.cursor()
        cur.execute("INSERT into users (username, type_id) " +
                    "VALUES (%s, %s) ", (username, userTypeId))
        self.conn.commit()

        return True

    # Function: GetUserIdFromUsername
    # Get a user id from a given username.
    #
    # Parameters:
    #     username - The username of the user to get an id for.
    #
    # Returns:
    #     The id of the user if successful, None on failure.
    def GetUserIdFromUsername(self, username):
        cur = self.conn.cursor()
        cur.execute("SELECT id FROM users WHERE username = %s", (username,))
        result = cur.fetchone()

        if not result:
            self.ReportError("Username not found.")
            return None

        (user_id,) = result
        return user_id

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
    #     True if the current user is a high user type than the user specified. False otherwise.
    def isCurrentUserRoleHigherThanUser(self, user_id):
        """
        Checks if the current user has higher privileges that the specified user_id.
        """
        cur = self.conn.cursor()
        cur.execute("SELECT type_id FROM users WHERE id = %s", (user_id,))
        result = cur.fetchone()

        if not result:
            self.ReportError("User not found.")
            return None

        (type_id,) = result
        if self.isRootUser() or (self.currentUser.getTypeId() < type_id):
            return True
        else:
            return False

    def CopyFastaToCopyDir(fasta_file, genome_id):

        if self.genomeCopyDir is None:
            self.ReportError("Need to set the genome storage directory.")
            return False
        if not os.path.isdir(self.genomeCopyDir):
            self.ReportError("Genome storage directory is not a directory.")
            return False

        cur.execute("SELECT genome_source_id, id_at_source, external_id_prefix " +
                    "FROM genomes, genome_sources " +
                    "WHERE id = %s "+
                    "AND genome_source_id = genome_sources.id", (genome_id,))

        result = cur.fetchall()
        if len(result) == 0:
            self.ReportError("Genome id not found: %s." % genome_id)
            return False

        for (genome_source_id, id_at_source, external_id_prefix) in result:
            target_file_name = external_id_prefix + "_" + str(id_at_source)
            try:
                shutil.copy(fasta_file, os.path.join(self.genomeCopyDir, target_file_name))
            except:
                self.ReportError("Copy to genome storage dir failed.")
                return False
            return target_file_name

    def GenerateTempTableName(self):

        rng = random.SystemRandom()
        suffix = ''
        for i in range(0,10):
            suffix += rng.choice('abcefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
        return "TEMP" + suffix + str(int(time.time()))


    #def AddFastaGenome(self, fasta_file, copy_fasta, name, desc, force_overwrite=False, source=None, id_at_source=None, completeness=0, contamination=0):
    #    cur = self.conn.cursor()
    #
    #    genome_id = self.AddFastaGenomeWorking(cur, fasta_file, name, desc, force_overwrite, source, id_at_source, completeness, contamination)
    #    if (genome_id):
    #        if copy_fasta:
    #            target_file_name = self.CopyFastaToCopyDir(fasta_file, genome_id)
    #            if target_file_name:
    #                self.conn.commit()
    #                return True
    #        else:
    #            self.conn.commit()
    #            return True
    #
    #    self.conn.rollback()
    #    return False

    def AddManyFastaGenomes(self, batchfile, checkM_file, modify_genome_list_id=None,
                            new_genome_list_name=None, force_overwrite=False):

        checkm_fh = open(checkM_file, "rb")

        expected_headers = ["Bin Id", "Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2",
                            "3", "4", "5+", "Completeness", "Contamination",  "Strain heterogeneity"]
        try:
        # Check the CheckM headers are consistent
            if checkm_fh.readline().rstrip() != "\t".join(expected_headers):
                raise GenomeDatabaseError("CheckM file header inconsistent. Expected: \n\t%s\nGot:\n\t%s\n" %
                                          ("\t".join(expected_headers), checkm_fh.readline().rstrip()))

            # Populate CheckM results dict
            checkM_results_dict = {}

            for line in checkm_fh:
                line = line.rstrip()
                splitline = line.split("\t")
                file_name, completeness, contamination = splitline[0], splitline[11], splitline[12]

                checkM_results_dict[file_name] = {"completeness" : completeness, "contamination" : contamination}

            checkm_fh.close()

            cur = self.conn.cursor()

            if modify_genome_list_id is not None:
                if new_genome_list_name is not None:
                    raise GenomeDatabaseError("Unable to both modify and create genome lists at the same time.")
                if self.GetGenomeIdListFromGenomeListId(modify_genome_list_id) is None:
                    raise GenomeDatabaseError("Unable to add genomes to list %s." % modify_genome_list_id)

            if new_genome_list_name is not None:
                owner_id = None
                if not self.currentUser.isRootUser():
                    owner_id = self.currentUser.getUserId()
                modify_genome_list_id = self.CreateGenomeListWorking(cur, [], new_genome_list_name, "", owner_id)
                if modify_genome_list_id is None:
                    raise GenomeDatabaseError("Unable to create the new genome list.")

            # Add the genomes
            fasta_paths_to_copy = []
            added_genome_ids = []

            fh = open(batchfile, "rb")
            for line in fh:
                line = line.rstrip()
                splitline = line.split("\t")
                if len(splitline) < 5:
                    splitline += [None] * (5 - len(splitline))
                (fasta_path, name, desc, source_name, id_at_source) = splitline

                abs_path = os.path.abspath(fasta_path)
                basename = os.path.splitext(os.path.basename(abs_path))[0]

                if basename not in checkM_results_dict:
                    raise GenomeDatabaseError("Couldn't find checkM result for %s (%s)" % (name,abs_path))

                genome_id = self.AddFastaGenomeWorking(
                    cur, abs_path, name, desc, modify_genome_list_id, force_overwrite, source_name, id_at_source,
                    checkM_results_dict[basename]["completeness"], checkM_results_dict[basename]["contamination"]
                )

                # Rollback everything if addition fails
                if not (genome_id):
                    raise GenomeDatabaseError("Failed to add genome: %s" % abs_path)

                added_genome_ids.append(genome_id)

                fasta_paths_to_copy.append(abs_path)

            if not self.EditGenomeListWorking(cur, modify_genome_list_id, genome_ids=added_genome_ids, operation='add'):
                raise GenomeDatabaseError("Unable to add genomes to genome list.")

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return False
        except:
            self.conn.rollback()
            raise

        self.conn.commit()
        return True

    # Function: AddFastaGenomeWorking
    # Checks if the current user is a higher user type than the specified user.
    #
    # Parameters:
    #     genome_list_id - None: don't add to genome list. Number: add to existing genome list id
    #
    # Returns:
    #     Genome id if it was added. None if it fails.
    def AddFastaGenomeWorking(self, cur, fasta_file_path, name, desc, genome_list_id=None, force_overwrite=False,
                              source=None, id_at_source=None, completeness=0, contamination=0):
        try:
            try:
                fasta_fh = open(fasta_file_path, "rb")
            except:
                raise GenomeDatabaseError("Cannot open Fasta file: " + fasta_file_path)

            m = hashlib.sha256()
            for line in fasta_fh:
                m.update(line)
            fasta_sha256_checksum = m.hexdigest()
            fasta_fh.close()

            if source is None:
                source = self.defaultGenomeSourceName


            if genome_list_id is not None:
                if self.GetGenomeIdListFromGenomeListId(genome_list_id) is None:
                    raise GenomeDatabaseError("Unable to add genome to list %s." % genome_list_id)

            cur.execute("SELECT id, external_id_prefix, user_accessible FROM genome_sources WHERE name = %s" , (source,))
            source_id = None
            prefix = None

            for (id, external_id_prefix, user_accessible) in cur:
                if (not user_accessible):
                    if id_at_source == None:
                        raise GenomeDatabaseError("Cannot auto generate ids at source for the %s genome source." % source)
                    if (not self.currentUser.isRootUser()):
                        raise GenomeDatabaseError("Only the root user can add genomes to the %s genome source." % source)
                source_id = id
                prefix = external_id_prefix
                break

            if source_id is None:
                raise GenomeDatabaseError("Could not find the %s genome source." % source)

            if id_at_source is None:
                cur.execute("SELECT id_at_source FROM genomes WHERE genome_source_id = %s order by id_at_source::int desc", (source_id,))
                last_id = None
                for (last_id_at_source, ) in cur:
                    last_id = last_id_at_source
                    break

                # Generate a new id (for user-accessible lists only)
                if (last_id is None):
                    new_id = 1
                else:
                    new_id = int(last_id) + 1

                if id_at_source is None:
                    id_at_source = str(new_id)

            added = datetime.datetime.now()

            owner_id = None
            if not self.currentUser.isRootUser():
                owner_id = self.currentUser.getUserId()

            cur.execute("SELECT id FROM genomes WHERE genome_source_id = %s AND id_at_source = %s", (source_id, id_at_source))

            result = cur.fetchall()

            columns = "(name, description, owned_by_root, owner_id, fasta_file_location, " + \
                      "fasta_file_sha256, genome_source_id, id_at_source, date_added, checkm_completeness, checkm_contamination)"


            if len(result):
                genome_id = result[0]
                if force_overwrite:
                    raise GenomeDatabaseError("Force overwrite not implemented yet")
                else:
                    raise GenomeDatabaseError("Genome source '%s' already contains id '%s'. Use -f to force an overwrite." % (source, id_at_source))
            else:
                cur.execute("INSERT INTO genomes " + columns + " "
                            "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) " +
                            "RETURNING id" ,
                            (name, desc, self.currentUser.isRootUser(), owner_id, fasta_file_path, fasta_sha256_checksum, source_id, id_at_source, added, completeness, contamination))


            genome_id = cur.fetchone()[0]

            # TODO: Add to genome list if required

            return genome_id

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None
        except:
            raise

    def ExternalGenomeIdsToGenomeIds(self, external_ids):
        try:
            cur = self.conn.cursor()

            map_sources_to_ids = {}

            for external_id in external_ids:
                try:
                    (source_prefix, id_at_source) = external_id.split("_", 1)
                except ValueError:
                    raise GenomeDatabaseError("All genome ids must have the form <prefix>_<id>. Offending id: %s" % str(external_id))

                if source_prefix not in map_sources_to_ids:
                    map_sources_to_ids[source_prefix] = {}
                map_sources_to_ids[source_prefix][id_at_source] = external_id

            temp_table_name = self.GenerateTempTableName()

            if len(map_sources_to_ids.keys()):
                cur.execute("CREATE TEMP TABLE %s (prefix text)" % (temp_table_name,) )
                query = "INSERT INTO {0} (prefix) VALUES (%s)".format(temp_table_name)
                cur.executemany(query, [(x,) for x in map_sources_to_ids.keys()])
            else:
                raise GenomeDatabaseError("No genome sources found for these ids. %s" % str(external_ids))

            # Find any given tree prefixes that arent in the genome sources
            query = ("SELECT prefix FROM {0} " +
                     "WHERE prefix NOT IN ( " +
                        "SELECT external_id_prefix " +
                        "FROM genome_sources)").format(temp_table_name)

            cur.execute(query)

            missing_genome_sources = {}
            for (query_prefix,) in cur:
                missing_genome_sources[query_prefix] = map_sources_to_ids[query_prefix].values()

            if len(missing_genome_sources.keys()):
                errors = []
                for (source_prefix, offending_ids) in missing_genome_sources.items():
                    errors.append("(%s) %s" % (source_prefix, str(offending_ids)))
                raise GenomeDatabaseError("Cannot find the relevant genome source id for the following ids, check the IDs are correct: " +
                                          ", ".join(errors))

            # All genome sources should be good, find ids
            result_ids = []
            for source_prefix in map_sources_to_ids.keys():

                # Create a table of requested external ids from this genome source
                temp_table_name = self.GenerateTempTableName()
                cur.execute("CREATE TEMP TABLE %s (id_at_source text)" % (temp_table_name,) )
                query = "INSERT INTO {0} (id_at_source) VALUES (%s)".format(temp_table_name)
                cur.executemany(query, [(x,) for x in map_sources_to_ids[source_prefix].keys()])

                # Check to see if there are any that don't exist
                query = ("SELECT id_at_source FROM {0} " +
                         "WHERE id_at_source NOT IN ( " +
                            "SELECT id_at_source " +
                            "FROM genomes, genome_sources " +
                            "WHERE genome_source_id = genome_sources.id "+
                            "AND external_id_prefix = %s)").format(temp_table_name)

                cur.execute(query, (source_prefix,))

                missing_ids = []
                for (id_at_source, ) in cur:
                    missing_ids.append(source_prefix + "_" + id_at_source)

                if missing_ids:
                    raise GenomeDatabaseError("Cannot find the the following genome ids, check the IDs are correct: %s" % str(missing_ids))

                # All exist, so get their ids.
                query = ("SELECT genomes.id FROM genomes, genome_sources " +
                         "WHERE genome_source_id = genome_sources.id "+
                         "AND id_at_source IN ( " +
                            "SELECT id_at_source " +
                            "FROM {0} )"+
                         "AND external_id_prefix = %s").format(temp_table_name)

                cur.execute(query, (source_prefix,))

                for (genome_id, ) in cur:
                    result_ids.append(genome_id)

            return result_ids

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None
        except:
            raise

    def GetAllGenomeIds(self):
        try:
            cur = self.conn.cursor()

            query = "SELECT id FROM genomes";
            cur.execute(query)

            result_ids = []
            for (genome_id, ) in cur:
                result_ids.append(genome_id)

            return result_ids

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None
        except:
            raise

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
                        raise GenomeDatabaseError("Cannot open batchfile: " + batchfile)

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
        except:
            raise

    def PrintGenomesDetails(self, genome_id_list):
        try:
            if not genome_id_list:
                raise GenomeDatabaseError("Unable to print genomes. No genomes found.")
            
            cur = self.conn.cursor()

            columns = "genomes.id, genomes.name, description, owned_by_root, username, fasta_file_location, " + \
                       "external_id_prefix || '_' || id_at_source as external_id, date_added, checkm_completeness, checkm_contamination"

            cur.execute("SELECT " + columns + " FROM genomes " +
                        "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
                        "JOIN genome_sources AS sources ON genome_source_id = sources.id " +
                        "AND genomes.id in %s "+
                        "ORDER BY genomes.id ASC", (tuple(genome_id_list),))

            print "\t".join(("genome_id", "name", "description", "owner", "fasta", "data_added", "completeness", "contamination"))

            for (genome_id, name, description, owned_by_root, username, fasta_file_location,
                 external_id, date_added, completeness, contamination) in cur:
                print "\t".join(
                    (external_id, name, description, ("(root)" if owned_by_root else username),
                     fasta_file_location, str(date_added), str(completeness), str(contamination))
                )
            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False
        except:
            raise

    def AddMarkers(self, batchfile, modify_marker_set_id=None, new_marker_set_name=None,
                   force_overwrite=False):

        # Add the markers
        cur = self.conn.cursor()
        marker_paths_to_copy = []

        fh = open(batchfile, "rb")
        for line in fh:
            line = line.rstrip()
            splitline = line.split("\t")
            if len(splitline) < 5:
                splitline += [None] * (5 - len(splitline))
            (marker_path, name, desc, database_name, database_specific_id) = splitline

            abs_path = os.path.abspath(marker_path)

            marker_id = self.AddMarkersWorking(cur, abs_path, name, desc, force_overwrite, database_name, database_specific_id)

            # Rollback everything if addition fails
            if not (marker_id):
                self.conn.rollback()
                return False

            marker_paths_to_copy.append(abs_path)

        #if copy_fastas:
        #    # TODO: Copy the fastas if required, rollback if fails
        #    pass

        self.conn.commit()
        return True

    def AddMarkersWorking(self, cur, marker_file_path, name, desc, marker_set_id=None, force_overwrite=False,
                          database=None, database_specific_id=None):
        try:
            marker_fh = open(marker_file_path, "rb")
        except:
            self.ReportError("Cannot open Marker file: " + marker_file_path)
            marker_fh.close()
            return None

        seen_name_line = False
        model_length = None

        m = hashlib.sha256()
        for line in marker_fh:
            if line[:4] == 'NAME':
                if seen_name_line:
                    self.ReportError("Marker file contains more than one model. Offending file: " + marker_file_path)
                    return None
                seen_name_line = True
            elif line[:4] == 'LENG':
                try:
                    model_length = int(line[4:])
                except:
                    self.ReportError("Unable to convert model length into integer value. Offending line: %s. Offending file %s." % (line, marker_file_path))
                    return None
            m.update(line)

        if model_length is None:
            self.ReportError("Model file does not give specify marker length. Offending file %s." % marker_file_path)
            return None

        if model_length <= 0:
            self.ReportError("Model file specifies invalid marker length. Length: %i. Offending file %s." % (model_length, marker_file_path))
            return None

        if marker_set_id is not None:
            if self.GetMarkerIdListFromMarkerSetId(marker_set_id) is None:
                raise GenomeDatabaseError("Unable to add marker to set %s." % marker_set_id)

        marker_sha256_checksum = m.hexdigest()
        marker_fh.close()

        if database is None:
            database = self.defaultMarkerDatabaseName

        cur.execute("SELECT id, user_accessible FROM marker_databases WHERE name = %s" , (database,))
        database_id = None

        for (id, user_accessible) in cur:
            if (not user_accessible):
                if database_specific_id == None:
                    self.ReportError("Cannot auto generate database specific ids for the %s marker database." % database)
                    return None
                if (not self.currentUser.isRootUser()):
                    self.ReportError("Only the root user can add markers to the %s database." % database)
                    return None
            database_id = id
            break

        if database_id is None:
            self.ReportError("Could not find the %s marker databasee." % database)
            return None

        if database_specific_id is None:

            cur.execute("SELECT database_specific_id FROM markers WHERE database_id = %s order by database_specific_id::int desc", (database_id,))
            last_id = None
            for (last_database_specific_id, ) in cur:
                last_id = last_database_specific_id
                break

            # Generate a new id (for user-accessible lists only)
            if (last_id is None):
                new_id = 1
            else:
                new_id = int(last_id) + 1

            if database_specific_id is None:
                database_specific_id = str(new_id)

        owner_id = None
        if not self.currentUser.isRootUser():
            owner_id = self.currentUser.getUserId()

        cur.execute("SELECT id FROM markers WHERE database_id = %s AND database_specific_id = %s", (database_id, database_specific_id))

        result = cur.fetchall()

        columns = "(name, owned_by_root, owner_id, marker_file_location, " + \
                  "marker_file_sha256, database_id, database_specific_id, size)"


        if len(result):
            marker_id = result[0]
            if force_overwrite:
                self.ReportError("Force overwrite not implemented yet")
                return None
            else:
                self.ReportError("Marker database '%s' already contains id '%s'. Use -f to force an overwrite." % (database, database_specific_id))
                return None
        else:
            cur.execute("INSERT INTO markers " + columns + " "
                        "VALUES (%s, %s, %s, %s, %s, %s, %s, %s) " +
                        "RETURNING id" ,
                        (name, self.currentUser.isRootUser(), owner_id, marker_file_path, marker_sha256_checksum, database_id, database_specific_id, model_length))

        marker_id = cur.fetchone()[0]

        # TODO: Add to marker set if needed

        return marker_id


    # Function: GetMarkerIdListFromMarkerSetId
    # Given a marker set id, return all the ids of the markers contained within that marker set.
    #
    # Parameters:
    #     marker_set_id - The marker set id of the marker set whose contents needs to be retrieved.
    #
    # Returns:
    #   A list of all the marker ids contained within the specified marker set, None on failure.
    def GetMarkerIdListFromMarkerSetId(self, marker_set_id):

        cur = self.conn.cursor()

        cur.execute("SELECT id, owner_id, owned_by_root, private " +
                    "FROM marker_sets " +
                    "WHERE id = %s ", (marker_set_id,))

        result = cur.fetchone()

        if not result:
            self.ReportError("No marker set with id: %s" % str(marker_set_id))
            return None
        else:
            (list_id, owner_id, owned_by_root, private) = result
            if private and (not self.currentUser.isRootUser()) and (owned_by_root or owner_id != self.currentUser.getUserId()):
                self.ReportError("Insufficient permission to view marker set: %s" % str(marker_set_id))
                return None


        cur.execute("SELECT marker_id " +
                    "FROM marker_set_contents " +
                    "WHERE set_id = %s ", (marker_set_id,))

        return [marker_id for (marker_id,) in cur.fetchall()]


    def FindUncalculatedMarkersForGenomeId(self, genome_id, marker_ids):
        
        cur = self.conn.cursor()
        
        cur.execute("SELECT marker_id, sequence " +
                    "FROM aligned_markers " +
                    "WHERE genome_id = %s ", (genome_id,))
        
        marker_id_dict = dict(cur.fetchall())
        
        return [x for x in marker_ids if x not in marker_id_dict]
                

    def MakeTreeData(self, marker_ids, genome_ids, directory, prefix, profile=None, config_dict=None, build_tree=True):

        cur = self.conn.cursor()

        if profile is None:
            profile = profiles.ReturnDefaultProfileName()
        if profile not in profiles.profiles:
            self.ReportError("Unknown Profile: " + profile)
            return None
        if not(os.path.exists(directory)):
            os.makedirs(directory)

        uncalculated_marker_dict = {}
        uncalculated_marker_count = 0

        for genome_id in genome_ids:
            uncalculated = self.FindUncalculatedMarkersForGenomeId(genome_id, marker_ids)
            if len(uncalculated) != 0:
                uncalculated_marker_dict[genome_id] = uncalculated
                uncalculated_marker_count += len(uncalculated)

        print "%i genomes contain %i uncalculated markers." % (len(uncalculated_marker_dict.keys()), uncalculated_marker_count)
        # TODO: Add Confirm step

        for (genome_id, uncalculated) in incalculated_marker_dict:
            self.RecalculateMarkersForGenome(genome_id, uncalculated)

        return profiles.profiles[profile].MakeTreeData(self, marker_ids, genome_ids,
                                                       directory, prefix, config_dict)

    # Function: CreateGenomeListWorking
    # Creates a new genome list in the database
    #
    # Parameters:
    #     cur -
    #     genome_id_list - A list of genome ids to add to the new list.
    #     name - The name of the newly created list.
    #     description - A description of the newly created list.
    #     owner_id - The id of the user who will own this list.
    #     private - Bool that denotes whether this list is public or private.
    #
    # Returns:
    #    The genome list id of the newly created list.
    def CreateGenomeListWorking(self, cur, genome_id_list, name, description, owner_id=None, private=True):

        if (owner_id is None):
            if not self.currentUser.isRootUser():
                self.ReportError("Only the root user can create root owned lists.")
                return None
        else:
            if (not self.currentUser.isRootUser()) and (self.currentUser.getUserId() != owner_id):
                self.ReportError("Only the root user may create lists on behalf of other people.")
                return None

        query = "INSERT INTO genome_lists (name, description, owned_by_root, owner_id, private) VALUES (%s, %s, %s, %s, %s) RETURNING id"
        cur.execute(query, (name, description, owner_id is None, owner_id, private))
        (genome_list_id, ) = cur.fetchone()

        query = "INSERT INTO genome_list_contents (list_id, genome_id) VALUES (%s, %s)"
        cur.executemany(query, [(genome_list_id, x) for x in genome_id_list])

        return genome_list_id


    # Function: GetGenomeIdListFromGenomeListId
    # Given a genome list id, return all the ids of the genomes contained within that genome list.
    #
    # Parameters:
    # genome_list_id - The genome list id of the genome list whose contents needs to be retrieved.
    #
    # Returns:
    # A list of all the genome ids contained within the specified genome list, None on failure.
    def GetGenomeIdListFromGenomeListId(self, genome_list_id):
        try:
            cur = self.conn.cursor()

            cur.execute("SELECT id, owner_id, owned_by_root, private " +
                        "FROM genome_lists " +
                        "WHERE id = %s ", (genome_list_id,))
            result = cur.fetchone()

            if not result:
                raise GenomeDatabaseError("No genome list with id: %s" % str(genome_list_id))
            else:
                (list_id, owner_id, owned_by_root, private) = result
                if private and (not self.currentUser.isRootUser()) and (owned_by_root or owner_id != self.currentUser.getUserId()):
                    raise GenomeDatabaseError("Insufficient permission to view genome list: %s" % str(genome_list_id))

            cur.execute("SELECT genome_id " +
                        "FROM genome_list_contents " +
                        "WHERE list_id = %s", (genome_list_id,))

            return [genome_id for (genome_id,) in cur.fetchall()]

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None
        except:
            raise

    # Function: GetGenomeIdListFromGenomeListIds
    # Given a list of ids, return all the ids of the genomes contained
    # within that genome list.
    #
    # Parameters:
    #     genome_list_ids - A list of genome list ids whose contents needs to be retrieved.
    #
    # Returns:
    #   A list of all the genome ids contained within the specified genome list(s), None on failure.
    def GetGenomeIdListFromGenomeListIds(self, genome_list_ids):
        try:
            cur = self.conn.cursor()

            temp_table_name = self.GenerateTempTableName()

            if genome_list_ids:
                cur.execute("CREATE TEMP TABLE %s (id integer)" % (temp_table_name,) )
                query = "INSERT INTO {0} (id) VALUES (%s)".format(temp_table_name)
                cur.executemany(query, [(x,) for x in genome_list_ids])
            else:
                raise GenomeDatabaseError("No genome lists given. Can not retrieve IDs" )

            # Find any ids that don't have genome lists
            query = ("SELECT id FROM {0} " +
                     "WHERE id NOT IN ( " +
                        "SELECT id " +
                        "FROM genome_lists)").format(temp_table_name)

            cur.execute(query)

            missing_list_ids = []
            for (list_id,) in cur:
                missing_list_ids.append(list_id)

            if missing_list_ids:
                raise GenomeDatabaseError("Unknown genome list id(s) given. %s" % str(missing_list_ids))

            # Find any genome list ids that we dont have permission to view
            cur.execute("SELECT id, owner_id, owned_by_root, private " +
                        "FROM genome_lists " +
                        "WHERE id in %s ", (tuple(genome_list_ids),))

            no_permission_list_ids = []
            for (list_id, owner_id, owned_by_root, private) in cur:
                if private and (not self.currentUser.isRootUser()) and (owned_by_root or owner_id != self.currentUser.getUserId()):
                    no_permission_list_ids.append(list_id)

            if no_permission_list_ids:
                raise GenomeDatabaseError("Insufficient permission to view genome lists: %s" % str(no_permission_list_ids))

            cur.execute("SELECT genome_id " +
                        "FROM genome_list_contents " +
                        "WHERE list_id in %s", (tuple(genome_list_ids),))

            return [genome_id for (genome_id,) in cur.fetchall()]

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None
        except:
            raise


    # Function: GetVisibleGenomeLists
    # Get all the genome lists that the current user can see.
    #
    # Parameters:
    #     owner_id - Get visible genome lists owned by this user with this id. If not specified, get all root owned lists.
    #
    # Returns:
    #   A list containing a tuple for each visible genome list. The tuple contains the genome list id, genome list name, genome list description,
    # and username of the owner of the list (id, name, description, username).
    def GetVisibleGenomeListsByOwner(self, owner_id=None):
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
            conditional_query += "AND (private = False OR owner_id = %s)"
            params.append(self.currentUser.getUserId())

        cur.execute("SELECT id " +
                    "FROM genome_lists " +
                    "WHERE 1 = 1 " +
                    conditional_query, params)

        return [list_id for (list_id,) in cur]

    def GetAllVisibleGenomeListIds(self):
        cur = self.conn.cursor()

        conditional_query = ""
        params = []

        if not self.currentUser.isRootUser():
            conditional_query += "AND (private = False OR owner_id = %s)"
            params.append(self.currentUser.getUserId())

        cur.execute("SELECT id " +
                    "FROM genome_lists " +
                    "WHERE 1 = 1 " +
                    conditional_query, params)

        return [list_id for (list_id,) in cur]

    def ViewGenomeListsContents(self, list_ids):
        try:
            genome_id_list = self.GetGenomeIdListFromGenomeListIds(list_ids)

            if genome_id_list is None:
                raise GenomeDatabaseError("Unable to view genome list. Can not retrieve genomes IDs for lists: %s" % str(list_ids))

            if not self.PrintGenomesDetails(genome_id_list):
                raise GenomeDatabaseError("Unable to view genome list. Printing to screen failed of genome ids. %s" % str(genome_id_list))

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False
        except:
            raise

    def PrintGenomeListsDetails(self, genome_list_ids):
        try:
            cur = self.conn.cursor()

            if not genome_list_ids:
                raise GenomeDatabaseError("Unable to print genome details: No genomes given." )
            
            if not self.currentUser.isRootUser():
                cur.execute("SELECT id " +
                            "FROM genome_lists as lists " +
                            "WHERE lists.private = True " +
                            "AND lists.id in %s " +
                            "AND (owned_by_root = True OR owner_id != %s)", (tuple(genome_list_ids), self.currentUser.getUserId()))

                unviewable_list_ids = [list_id for (list_id, ) in cur]
                if unviewable_list_ids:
                    raise GenomeDatabaseError("Insufficient privileges to view genome lists: %s." % str(unviewable_list_ids))


            cur.execute(
                "SELECT lists.id, lists.name, lists.description, lists.private, lists.owned_by_root, users.username, count(contents.list_id) " +
                "FROM genome_lists as lists " +
                "LEFT OUTER JOIN users ON lists.owner_id = users.id " +
                "JOIN genome_list_contents as contents ON contents.list_id = lists.id " +
                "WHERE lists.id in %s " +
                "GROUP by lists.id, users.username " +
                "ORDER by lists.id asc " , (tuple(genome_list_ids),)
            )

            print "\t".join(("list_id", "name", "description", "owner", "visibility", "genome_count"))

            for (list_id, name, description, private, owned_by_root, username, genome_count) in cur:
                print "\t".join(
                    (str(list_id), name, description, ("(root)" if owned_by_root else username), ("private" if private else "public"), str(genome_count))
                )
            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False
        except:
            raise

    def HasPermissionToViewGenomeList(self, genome_list_id):
        try:
            cur = self.conn.cursor()
        
            cur.execute("SELECT owner_id, owned_by_root, private " +
                        "FROM genome_lists " +
                        "WHERE id = %s ", (genome_list_id,))
        
            result = cur.fetchone()
            
            if not result:
                raise GenomeDatabaseError("No genome list with id: %s" % str(genome_list_id))
            
            (owner_id, owned_by_root) = result
            
            if not self.currentUser.isRootUser():
                if private and (owned_by_root or owner_id != self.currentUser.getUserId()):
                    return False
            else:
                if not owned_by_root:
                    self.ReportError("Root user editing of other users lists not yet implmented.")
                    return False
            
            return True
            
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return None
        except:
            raise
        
        if not self.currentUser.isRootUser():
            conditional_query += "AND (private = False OR owner_id = %s)"
            params.append(self.currentUser.getUserId())

        cur.execute("SELECT id " +
                    "FROM genome_lists " +
                    "WHERE 1 = 1 " +
                    conditional_query, params)

    def HasPermissionToEditGenomeList(self, genome_list_id):
        try:
            cur = self.conn.cursor()
        
            cur.execute("SELECT owner_id, owned_by_root " +
                        "FROM genome_lists " +
                        "WHERE id = %s ", (genome_list_id,))
        
            result = cur.fetchone()
            
            if not result:
                raise GenomeDatabaseError("No genome list with id: %s" % str(genome_list_id))
            
            (owner_id, owned_by_root) = result
            
            if not self.currentUser.isRootUser():
                if owned_by_root or owner_id != self.currentUser.getUserId():
                    return False
            else:
                if not owned_by_root:
                    self.ReportError("Root user editing of other users lists not yet implmented.")
                    return False
            
            return True
            
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return None
        except:
            raise
        
        if not self.currentUser.isRootUser():
            conditional_query += "AND (private = False OR owner_id = %s)"
            params.append(self.currentUser.getUserId())

        cur.execute("SELECT id " +
                    "FROM genome_lists " +
                    "WHERE 1 = 1 " +
                    conditional_query, params)
    
    def EditGenomeList(self, genome_list_id, batchfile=None, genomes_external_ids=None, operation=None, name=None, description=None, private=None):        
        
        cur = self.conn.cursor()
    
        if batchfile:
            if genomes_external_ids is None:
                genomes_external_ids = []
            for line in fh:
                line = line.rstrip()
                genomes_external_ids.append(line)
                
        if genomes_external_ids is not None:
            genomes_external_ids = self.ExternalGenomeIdsToGenomeIds(genomes_external_ids)

        if not self.EditGenomeListWorking(cur, genome_list_id, genomes_external_ids, operation, name, description, private):
            self.conn.rollback()
            return False
        
        self.conn.commit()
        return True
        
    
    def EditGenomeListWorking(self, cur, genome_list_id, genome_ids=None, operation=None, name=None, description=None, private=None):
        try:
            edit_permission = self.HasPermissionToEditGenomeList(genome_list_id)
            if edit_permission is None:
                raise GenomeDatabaseError("Unable to retrieve genome list id for editing. Offending list id: %s" % genome_list_id)
            elif edit_permission == False:
                raise GenomeDatabaseError("Insufficent permissions to edit this genome list. Offending list id: %s" % genome_list_id)
            
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
                cur.execute("UPDATE genome_lists SET " + update_query + " WHERE id = %s", params + [genome_list_id])
            
            temp_table_name = self.GenerateTempTableName()


            if operation is not None:
                
                if len(genome_ids) == 0:
                    raise GenomeDatabaseError("No genome ids given to perform '%s' operation." % operation)
                
                cur.execute("CREATE TEMP TABLE %s (id integer)" % (temp_table_name,) )
                query = "INSERT INTO {0} (id) VALUES (%s)".format(temp_table_name)
                cur.executemany(query, [(x,) for x in genome_ids])
        
                if operation == 'add':
                    query = ("INSERT INTO genome_list_contents (list_id, genome_id) " +
                             "SELECT %s, id FROM {0} " +
                             "WHERE id NOT IN ( " +
                                "SELECT genome_id " +
                                "FROM genome_list_contents " +
                                "WHERE list_id = %s)").format(temp_table_name)
                    cur.execute(query, (genome_list_id, genome_list_id))
                elif operation == 'remove':
                    query = ("DELETE FROM genome_list_contents " +
                            "WHERE list_id = %s " +
                            "AND genome_id IN ( " +
                                "SELECT id " +
                                "FROM {0})").format(temp_table_name)
                    cur.execute(query, [genome_list_id])
                else:
                    raise GenomeDatabaseError("Unknown genome list edit operation: %s" % operation)
            
            return True
            
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None
        except:
            raise
    
