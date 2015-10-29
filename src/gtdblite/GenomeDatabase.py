import shutil
import os
import stat
import datetime
import hashlib
import logging
from multiprocessing import Pool

import prettytable

from gtdblite import Config
from gtdblite import ConfigMetadata
from gtdblite.User import User
from gtdblite.GenomeDatabaseConnection import GenomeDatabaseConnection
from gtdblite import MarkerCalculation
from gtdblite.MetadataManager import MetadataManager
from gtdblite import Tools
from gtdblite.GenomeFilter import GenomeFilter
from gtdblite.Exceptions import GenomeDatabaseError

from biolib.checksum import sha256


class GenomeDatabase(object):

    def __init__(self, threads=1):
        self.logger = logging.getLogger()

        self.conn = GenomeDatabaseConnection()
        self.currentUser = None
        self.errorMessages = []
        self.warningMessages = []
        self.debugMode = False

        self.threads = threads
        self.pool = Pool(threads)

        self.genomeFileSuffix = "_genomic.fna"
        self.proteinFileSuffix = "_protein.faa"
        self.ntGeneFileSuffix = "_protein.fna"
        self.gffFileSuffix = "_protein.gff"
        self.tranTableFileSuffix = "_translation_table.tsv"
        self.checksumSuffix = ".sha256"

        self.genomeCopyDir = None
        if Config.GTDB_GENOME_COPY_DIR:
            self.genomeCopyDir = Config.GTDB_GENOME_COPY_DIR

        self.markerCopyDir = None
        if Config.GTDB_MARKER_COPY_DIR:
            self.markerCopyDir = Config.GTDB_MARKER_COPY_DIR

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

    def CopyFastaToCopyDir(self, fasta_file, genome_id):
        try:
            if self.genomeCopyDir is None:
                raise GenomeDatabaseError(
                    "Need to set the genome storage directory.")

            if not os.path.isdir(self.genomeCopyDir):
                raise GenomeDatabaseError(
                    "Genome storage directory is not a directory.")

            cur.execute("SELECT genome_source_id, id_at_source, external_id_prefix " +
                        "FROM genomes, genome_sources " +
                        "WHERE id = %s " +
                        "AND genome_source_id = genome_sources.id", (genome_id,))

            result = cur.fetchall()
            if len(result) == 0:
                raise GenomeDatabaseError(
                    "Genome id not found: %s." % genome_id)

            (_genome_source_id, id_at_source, external_id_prefix) = result[0]

            target_file_name = external_id_prefix + "_" + str(id_at_source)
            target_file_path = os.path.join(
                self.genomeCopyDir, target_file_name)
            try:
                shutil.copy(fasta_file, target_file_path)
                os.chmod(
                    target_file_path, stat.S_IROTH | stat.S_IRGRP | stat.S_IRUSR)
            except:
                raise GenomeDatabaseError("Copy to genome storage dir failed.")

            return target_file_name

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # True on success, False otherwise (and on error)
    def AddManyFastaGenomes(self, batchfile, checkm_file, modify_genome_list_id=None,
                            new_genome_list_name=None, force_overwrite=False):
        try:
            self.logger.info("Reading CheckM file.")
            try:
                checkm_fh = open(checkm_file, "rb")
            except:
                raise GenomeDatabaseError(
                    "Cannot open CheckM file: " + checkm_file)

            required_headers = Tools.populate_required_headers(checkm_fh)

            # populate CheckM results dict
            checkm_results_dict = {}

            for line in checkm_fh:
                line = line.rstrip()
                splitline = line.split("\t")
                bin_id, completeness, contamination = (splitline[required_headers["Bin Id"]],
                                                       splitline[
                    required_headers["Completeness"]],
                    splitline[
                    required_headers["Contamination"]])
                lineage = None
                if required_headers.get("Marker lineage") is not None:
                    lineage = splitline[required_headers["Marker lineage"]]
                genome_count = None
                if required_headers.get("# genomes") is not None:
                    genome_count = splitline[required_headers["# genomes"]]
                marker_count = None
                if required_headers.get("# markers") is not None:
                    marker_count = splitline[required_headers["# markers"]]
                set_count = None
                if required_headers.get("# marker sets") is not None:
                    set_count = splitline[required_headers["# marker sets"]]
                heterogeneity = None
                if required_headers.get("Strain heterogeneity") is not None:
                    heterogeneity = splitline[
                        required_headers["Strain heterogeneity"]]

                checkm_results_dict[bin_id] = {"completeness": completeness,
                                               "contamination": contamination,
                                               "lineage": lineage,
                                               "genome_count": genome_count,
                                               "marker_count": marker_count,
                                               "set_count": set_count,
                                               "heterogeneity": heterogeneity}

            checkm_fh.close()

            cur = self.conn.cursor()

            if modify_genome_list_id is not None:
                if new_genome_list_name is not None:
                    raise GenomeDatabaseError(
                        "Unable to both modify and create genome lists at the same time.")

                has_permission = self.HasPermissionToEditGenomeList(
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

                modify_genome_list_id = self.CreateGenomeListWorking(
                    cur, [], new_genome_list_name, "", owner_id)
                if modify_genome_list_id is None:
                    raise GenomeDatabaseError(
                        "Unable to create the new genome list.")

            # add the genomes
            added_genome_dict = self.add_genome_list(
                cur, checkm_results_dict, batchfile, force_overwrite)

            # run Prodigal on Genomes having only a genome file
            self.logger.info(
                "Running Prodigal on genomes without identified genes.")
            fasta_paths_to_copy = Tools.runMultiProdigal(
                self.threads, added_genome_dict)

            if modify_genome_list_id is not None:
                if not self.EditGenomeListWorking(cur, modify_genome_list_id, genome_ids=added_genome_dict.keys(), operation='add'):
                    raise GenomeDatabaseError(
                        "Unable to add genomes to genome list.")

            copied_fasta_paths = []
            copied_genes_fasta_paths = []
            external_id_dict = {}

            cur.execute("SELECT genomes.id,user_editable, external_id_prefix || '_' || id_at_source as external_id " +
                        "FROM genomes, genome_sources " +
                        "WHERE genome_source_id = genome_sources.id " +
                        "AND genomes.id in %s", (tuple(added_genome_dict.keys()),))

            for (genome_id, user_editable, external_id) in cur:
                if user_editable:
                    external_id_dict[genome_id] = external_id

            if len(external_id_dict.keys()) > 0:
                username = None
                if self.currentUser.isRootUser():
                    username = self.currentUser.getElevatedFromUsername()
                else:
                    username = self.currentUser.getUsername()

                if username is None:
                    raise GenomeDatabaseError(
                        "Unable to determine user to add genomes under.")

                target_dir = os.path.join(self.genomeCopyDir, username)
                if os.path.exists(target_dir):
                    if not os.path.isdir(target_dir):
                        raise GenomeDatabaseError(
                            "Genome copy directory exists, but isn't a directory: %s" % (target_dir,))
                else:
                    os.mkdir(target_dir)

                try:
                    self.logger.info(
                        "Calculating and storing metadata for each genomes.")
                    metadata_mngr = MetadataManager()
                    for genome_id, external_id in external_id_dict.items():
                        sub_target_dir = os.path.join(target_dir, external_id)
                        if os.path.exists(sub_target_dir):
                            if not os.path.isdir(target_dir):
                                raise GenomeDatabaseError(
                                    "Genome copy directory exists, but isn't a directory: %s" % (target_dir,))
                        else:
                            os.mkdir(sub_target_dir)

                        target_file = os.path.join(
                            sub_target_dir, external_id + self.genomeFileSuffix)
                        shutil.copy(
                            fasta_paths_to_copy.get(genome_id).get('fasta_path'), target_file)
                        # os.chmod(target_file, stat.S_IROTH | stat.S_IRGRP | stat.S_IRUSR)
                        copied_fasta_paths.append(target_file)

                        genes_target_file = os.path.join(
                            sub_target_dir, external_id + self.proteinFileSuffix)
                        shutil.copy(
                            fasta_paths_to_copy.get(genome_id).get('aa_gene_file'), genes_target_file)
                        copied_genes_fasta_paths.append(genes_target_file)

                        checksum = sha256(genes_target_file)
                        fout = open(
                            genes_target_file + self.checksumSuffix, 'w')
                        fout.write(checksum)
                        fout.close()

                        nt_target_file = os.path.join(
                            sub_target_dir, external_id + self.ntGeneFileSuffix)
                        shutil.copy(
                            fasta_paths_to_copy.get(genome_id).get('nt_gene_file'), nt_target_file)

                        gff_target_file = os.path.join(
                            sub_target_dir, external_id + self.gffFileSuffix)
                        shutil.copy(
                            fasta_paths_to_copy.get(genome_id).get('gff_file'), gff_target_file)

                        trans_table_target_file = os.path.join(
                            sub_target_dir, external_id + self.tranTableFileSuffix)
                        shutil.copy(fasta_paths_to_copy.get(genome_id).get(
                            'translation_table_file'), trans_table_target_file)

                        cur.execute("UPDATE genomes SET fasta_file_location = %s , genes_file_location = %s WHERE id = %s", (
                            target_file, genes_target_file, genome_id))

                        # calculate metadata and store results in database
                        metadata_mngr.calculateMetadata(
                            target_file, gff_target_file, sub_target_dir)
                        metadata_mngr.storeMetadata(
                            sub_target_dir, genome_id, cur)

                except Exception as e:
                    try:
                        for copied_path in copied_fasta_paths:
                            os.unlink(copied_path)
                    except:
                        self.ReportWarning(
                            "Cleaning temporary copied files failed. May have orphan fastas in the genome copy directory.")
                    raise

            self.conn.commit()

            self.logger.info('Done.')

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return False
        except:
            self.conn.rollback()
            raise

    # Function: AddFastaGenomeWorking
    # Checks if the current user is a higher user type than the specified user.
    #
    # Parameters:
    #     genome_list_id - None: don't add to genome list. Number: add to existing genome list id
    #
    # Returns:
    #     Genome id if it was added. False if it fails.
    def AddFastaGenomeWorking(self, cur, fasta_file_path, name, desc, genome_list_id=None, force_overwrite=False,
                              source=None, id_at_source=None, gene_path=None, checkm_dict=None):
        try:
            fasta_sha256_checksum = sha256(fasta_file_path)

            gene_sha256_checksum = None
            if gene_path is not None:
                gene_sha256_checksum = sha256(gene_path)
            if source is None:
                source = self.defaultGenomeSourceName

            if genome_list_id is not None:
                has_permission = self.HasPermissionToEditGenomeList(
                    genome_list_id)
                if has_permission is None:
                    raise GenomeDatabaseError(
                        "Unable to add genome to list %s." % genome_list_id)
                elif not has_permission:
                    raise GenomeDatabaseError(
                        "Insufficient permission to add genome to genome list %s." % genome_list_id)

            cur.execute(
                "SELECT id, external_id_prefix, user_editable FROM genome_sources WHERE name = %s", (source,))
            source_id = None

            for (db_id, _external_id_prefix, user_editable) in cur:
                if (not user_editable):
                    if id_at_source is None:
                        raise GenomeDatabaseError(
                            "Cannot auto generate ids at source for the %s genome source." % source)
                    if (not self.currentUser.isRootUser()):
                        raise GenomeDatabaseError(
                            "Only the root user can add genomes to the %s genome source." % source)
                source_id = db_id
                break

            if source_id is None:
                raise GenomeDatabaseError(
                    "Could not find the %s genome source." % source)

            if id_at_source is None:
                cur.execute(
                    "SELECT id_at_source FROM genomes WHERE genome_source_id = %s order by id_at_source::int desc", (source_id,))
                last_id = None
                for (last_id_at_source,) in cur:
                    last_id = last_id_at_source
                    break

                cur.execute(
                    "SELECT last_auto_id FROM genome_sources WHERE id = %s ", (source_id,))
                for (last_auto_id,) in cur:
                    if last_id is None:
                        last_id = last_auto_id
                    else:
                        last_id = max(last_id, last_auto_id)
                    break

                # Generate a new id (for user editable lists only)
                if (last_id is None):
                    new_id = 1
                else:
                    new_id = int(last_id) + 1

                if id_at_source is None:
                    id_at_source = str(new_id)

                cur.execute(
                    "UPDATE genome_sources set last_auto_id = %s where id = %s", (new_id, source_id))

            added = datetime.datetime.now()

            owner_id = None
            if not self.currentUser.isRootUser():
                owner_id = self.currentUser.getUserId()

            cur.execute(
                "SELECT id FROM genomes WHERE genome_source_id = %s AND id_at_source = %s", (source_id, id_at_source))

            result = cur.fetchall()

            columns = "(name, description, owned_by_root, owner_id, fasta_file_location, " + \
                "fasta_file_sha256, genes_file_location, genes_file_sha256,genome_source_id, id_at_source, date_added)"

            if len(result):
                if force_overwrite:
                    raise GenomeDatabaseError(
                        "Force overwrite not implemented yet")
                else:
                    raise GenomeDatabaseError(
                        "Genome source '%s' already contains id '%s'. Use -f to force an overwrite." % (source, id_at_source))

            cur.execute("INSERT INTO genomes " + columns + " "
                        "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) " +
                        "RETURNING id",
                        (name, desc, self.currentUser.isRootUser(), owner_id, fasta_file_path, fasta_sha256_checksum, gene_path, gene_sha256_checksum, source_id, id_at_source, added))

            (genome_id,) = cur.fetchone()

            # We insert checkm information into the metadata_genes table
            column_checkm = "(id,checkm_completeness,checkm_contamination,checkm_marker_count,checkm_marker_lineage,checkm_genome_count,checkm_marker_set_count,checkm_strain_heterogeneity)"

            cur.execute("INSERT INTO metadata_genes " + column_checkm + " "
                        "VALUES (%s, %s, %s, %s, %s, %s, %s,%s) ", (genome_id, checkm_dict.get("completeness"),
                                                                    checkm_dict.get(
                            "contamination"),
                            checkm_dict.get(
                            "marker_count"),
                            checkm_dict.get(
                            "lineage"),
                            checkm_dict.get(
                            "genome_count"),
                            checkm_dict.get(
                            "set_count"),
                            checkm_dict.get("heterogeneity")))

            if genome_list_id:
                has_permission = self.HasPermissionToEditGenomeList(
                    genome_list_id)
                if has_permission is None:
                    raise GenomeDatabaseError(
                        "Error evaluating permission of genome list: %s", (genome_list_id,))
                elif not has_permission:
                    raise GenomeDatabaseError(
                        "Insufficient permission to edit genome list: %s", (genome_list_id,))
                cur.execute(
                    "INSERT INTO genome_list_contents (list_id, genome_id) VALUES (%s, %s)", (genome_list_id, genome_id))

            return genome_id

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def list_markers(self, cur=None, library=None):
        cur.execute("SELECT id_in_database FROM markers m " +
                    "LEFT OUTER JOIN marker_databases md ON m.marker_database_id = md.id " +
                    "WHERE md.name like '{0}'".format(library))
        listmarkers = [hmm_id for (hmm_id,) in cur.fetchall()]
        return listmarkers

    def add_genome_list(self, cur=None, checkm_results_dict=None, batchfile=None, force_overwrite=False):
        # Add the genomes
        added_genome_ids = {}
        fh = open(batchfile, "rb")
        for line in fh:
            line = line.rstrip()
            if line == '':
                self.ReportWarning(
                    "Encountered blank line in batch file. It has been ignored.")
                continue

            splitline = line.split("\t")
            if len(splitline) < 6:
                splitline += [None] * (6 - len(splitline))

            (fasta_path, name, desc, gene_path,
             source_name, id_at_source) = splitline

            if fasta_path is None or fasta_path == '':
                raise GenomeDatabaseError(
                    "Each line in the batch file must specify a path to the genome's fasta file.")

            if name is None or name == '':
                raise GenomeDatabaseError(
                    "Each line in the batch file must specify a name for the genome.")

            abs_path = os.path.abspath(fasta_path)
            bin_id = os.path.splitext(os.path.basename(abs_path))[0]

            abs_gene_path = None
            if gene_path is not None and gene_path != '':
                abs_gene_path = os.path.abspath(gene_path)

            if bin_id not in checkm_results_dict:
                raise GenomeDatabaseError(
                    "Couldn't find CheckM result for %s (%s). Bin identifier is %s." % (name, abs_path, bin_id))

            genome_id = self.AddFastaGenomeWorking(
                cur, abs_path, name, desc, None, force_overwrite, source_name, id_at_source, abs_gene_path,
                checkm_results_dict[bin_id])

            if not (genome_id):
                raise GenomeDatabaseError(
                    "Failed to add genome: %s" % abs_path)

            added_genome_ids[genome_id] = {
                "gene_path": abs_gene_path, "fasta_path": abs_path}
        return added_genome_ids

    # True if has permission. False if doesn't. None on error.
    def HasPermissionToEditGenome(self, genome_id):
        try:
            cur = self.conn.cursor()

            cur.execute("SELECT owner_id, owned_by_root "
                        "FROM genomes " +
                        "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
                        "AND genomes.id = %s", (genome_id,))

            result = cur.fetchone()

            if not result:
                raise GenomeDatabaseError(
                    "No genome list with id: %s" % genome_id)

            (owner_id, owned_by_root) = result

            if not self.currentUser.isRootUser():
                if (owned_by_root or owner_id != self.currentUser.getUserId()):
                    return False
            else:
                if not owned_by_root:
                    raise GenomeDatabaseError(
                        "Root user editing of other users genomes not yet implemented.")

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None

    # True if has permission. False if doesn't. None on error.
    def HasPermissionToEditGenomes(self, genome_ids):
        try:
            cur = self.conn.cursor()

            if not genome_ids:
                raise GenomeDatabaseError(
                    "Unable to retrieve genome permissions, no genomes given: %s" % str(genome_ids))

            cur.execute("SELECT genomes.id, owner_id, owned_by_root "
                        "FROM genomes " +
                        "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
                        "WHERE genomes.id in %s", (tuple(genome_ids),))

            for (genome_id, owner_id, owned_by_root) in cur:
                if not self.currentUser.isRootUser():
                    if (owned_by_root or owner_id != self.currentUser.getUserId()):
                        self.ReportWarning(
                            "Insufficient permissions to edit genome %s." % str(genome_id))
                        return False
                else:
                    if not owned_by_root:
                        self.ReportWarning(
                            "Root user editing of other users genomes not yet implemented.")
                        return False

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None

    # True on success. False on error/failure.
    def DeleteGenomes(self, batchfile=None, external_ids=None):
        '''
        Delete Genomes from either of the following parameters
        :param batchfile:test
        :param external_ids:

        '''
        try:
            cur = self.conn.cursor()

            if external_ids is None:
                external_ids = []

            if batchfile:
                fh = open(batchfile, "rb")
                for line in fh:
                    line = line.rstrip()
                    external_ids.append(line)

            genome_ids = self.ExternalGenomeIdsToGenomeIds(external_ids)

            if genome_ids is False:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Unable to retrieve genome ids.")

            has_permission = self.HasPermissionToEditGenomes(genome_ids)

            if has_permission is None:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Unable to retrieve permissions for genomes.")

            if has_permission is False:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Insufficient permissions.")

            if not self.Confirm("Are you sure you want to delete %i genomes (this action cannot be undone)" % len(genome_ids)):
                raise GenomeDatabaseError("User aborted database action.")

            paths_to_delete = []

            if self.genomeCopyDir is not None:
                cur.execute("SELECT fasta_file_location " +
                            "FROM genomes " +
                            "WHERE id in %s", (tuple(genome_ids),))

                for (fasta_path,) in cur:
                    # Check if path is a subdir of the copy dir
                    abs_dir = os.path.abspath(self.genomeCopyDir)
                    abs_file = os.path.abspath(fasta_path)

                    if abs_file.startswith(abs_dir):
                        paths_to_delete.append(fasta_path)

            cur.execute("DELETE FROM aligned_markers " +
                        "WHERE genome_id in %s", (tuple(genome_ids),))

            cur.execute("DELETE FROM genome_list_contents " +
                        "WHERE genome_id in %s", (tuple(genome_ids),))

            cur.execute("DELETE FROM genomes " +
                        "WHERE id in %s", (tuple(genome_ids),))

            try:
                for path_to_delete in paths_to_delete:
                    os.unlink(path_to_delete)
            except Exception as e:
                self.ReportWarning(
                    "Exception was raised when deleting genomes. Some orphans may remain. Exception message: %s" % e.message)

            self.conn.commit()
            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return False

    # List of genome ids on success. False on error.
    def ExternalGenomeIdsToGenomeIds(self, external_ids):
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
            self.ReportError(e.message)
            return False

    # List of genome ids on success. False on error.
    def GetAllGenomeIds(self):
        try:
            cur = self.conn.cursor()

            query = "SELECT id FROM genomes"
            cur.execute(query)

            result_ids = []
            for (genome_id,) in cur:
                result_ids.append(genome_id)

            return result_ids

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

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

    # True if success. False on failure/error.
    def PrintGenomesDetails(self, genome_id_list):
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
            table = prettytable.PrettyTable(("genome_id", "name", "description", "owner", "data_added"))
            table.align = 'l'
            table.hrules = prettytable.FRAME
            table.vrules = prettytable.NONE

            for (_genome_id, name, description, owned_by_root, username, external_id, date_added) in cur:
                tup = map(str, (external_id, name, description, ("root" if owned_by_root else username), date_added))
                table.add_row(tup)
            print table.get_string()

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # True if success. False on failure/error.
    def AddMarkers(self, batchfile, modify_marker_set_id=None, new_marker_set_name=None,
                   force_overwrite=False):
        try:
            cur = self.conn.cursor()

            if modify_marker_set_id is not None:
                if new_marker_set_name is not None:
                    raise GenomeDatabaseError(
                        "Unable to both modify and create marker sets at the same time.")
                has_permission = self.HasPermissionToEditMarkerSet(
                    modify_marker_set_id)
                if has_permission is None:
                    raise GenomeDatabaseError(
                        "Unable to add markers to set %s." % modify_marker_set_id)
                elif not has_permission:
                    raise GenomeDatabaseError(
                        "Insufficient permissions to add markers to set %s." % modify_marker_set_id)

            if new_marker_set_name is not None:
                owner_id = None
                if not self.currentUser.isRootUser():
                    owner_id = self.currentUser.getUserId()
                modify_marker_set_id = self.CreateMarkerSetWorking(
                    cur, [], new_marker_set_name, "", owner_id)
                if modify_marker_set_id is False:
                    raise GenomeDatabaseError(
                        "Unable to create the new marker set.")

            added_marker_ids = []

            fh = open(batchfile, "rb")
            for line in fh:
                line = line.rstrip()
                splitline = line.split("\t")
                if len(splitline) < 5:
                    splitline += [None] * (5 - len(splitline))
                (marker_path, name, desc, database_name,
                 id_in_database) = splitline

                abs_path = os.path.abspath(marker_path)

                marker_id = self.AddMarkerWorking(cur, abs_path, name, desc, None,
                                                  force_overwrite, database_name, id_in_database)

                if not (marker_id):
                    raise GenomeDatabaseError(
                        "Failed to add marker: %s" % abs_path)

                added_marker_ids.append(marker_id)

            if not self.EditMarkerSetWorking(cur, modify_marker_set_id, marker_ids=added_marker_ids, operation='add'):
                raise GenomeDatabaseError(
                    "Unable to add markers to marker set.")

            copied_hmm_paths = []
            hmm_paths_to_copy = {}

            cur.execute("SELECT markers.id, marker_file_location, user_editable, external_id_prefix || '_' || id_in_database as external_id " +
                        "FROM markers, marker_databases " +
                        "WHERE marker_database_id = marker_databases.id " +
                        "AND markers.id in %s", (tuple(added_marker_ids),))

            for (marker_id, abs_path, user_editable, external_id) in cur:
                if user_editable:
                    hmm_paths_to_copy[marker_id] = {'src_path': abs_path,
                                                    'external_id': external_id}

            if len(hmm_paths_to_copy.keys()) > 0:
                username = None
                if self.currentUser.isRootUser():
                    username = self.currentUser.getElevatedFromUsername()
                else:
                    username = self.currentUser.getUsername()

                if not username:
                    raise GenomeDatabaseError(
                        "Unable to determine user to add markers under.")

                target_dir = os.path.join(self.markerCopyDir, username)
                if os.path.exists(target_dir):
                    if not os.path.isdir(target_dir):
                        raise GenomeDatabaseError(
                            "Marker copy directory exists, but isn't a directory: %s" % (target_dir,))
                else:
                    os.mkdir(target_dir)

                try:
                    for (marker_id, details) in hmm_paths_to_copy.items():
                        target_file = os.path.join(
                            target_dir, details['external_id'] + ".hmm")
                        shutil.copy(details['src_path'], target_file)
                        copied_hmm_paths.append(target_file)

                        cur.execute(
                            "UPDATE markers SET marker_file_location = %s WHERE id = %s", (target_file, marker_id))

                except Exception as e:
                    try:
                        for copied_path in copied_hmm_paths:
                            os.unlink(copied_path)
                    except:
                        self.ReportWarning(
                            "Cleaning temporary copied files failed. May have orphan hmms in the marker copy directory.")
                    raise

            self.conn.commit()
            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return False
        except:
            self.conn.rollback()
            raise

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
                        last_id = max(last_id, last_auto_id)
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
            return marker_id

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def DeleteMarkers(self, batchfile=None, external_ids=None):
        try:
            cur = self.conn.cursor()

            if external_ids is None:
                external_ids = []

            if batchfile:
                fh = open(batchfile, "rb")
                for line in fh:
                    line = line.rstrip()
                    external_ids.append(line)

            marker_ids = self.ExternalMarkerIdsToMarkerIds(external_ids)
            if marker_ids is False:
                raise GenomeDatabaseError(
                    "Unable to delete markers. Unable to find markers.")

            has_permission = self.HasPermissionToEditMarkers(marker_ids)

            if has_permission is None:
                raise GenomeDatabaseError(
                    "Unable to delete markers. Unable to retrieve permissions for markers.")

            if has_permission is False:
                raise GenomeDatabaseError(
                    "Unable to delete markers. Insufficient permissions.")

            if not self.Confirm("Are you sure you want to delete %i markers (this action cannot be undone)" % len(marker_ids)):
                raise GenomeDatabaseError("User aborted database action.")

            paths_to_delete = []

            if self.genomeCopyDir is not None:
                cur.execute("SELECT marker_file_location " +
                            "FROM markers " +
                            "WHERE id in %s", (tuple(marker_ids),))

                for (hmm_path,) in cur:
                    # Check if path is a subdir of the copy dir
                    abs_dir = os.path.abspath(self.markerCopyDir)
                    abs_file = os.path.abspath(hmm_path)

                    if abs_file.startswith(abs_dir):
                        paths_to_delete.append(hmm_path)

            cur.execute("DELETE FROM marker_set_contents " +
                        "WHERE marker_id in %s", (tuple(marker_ids),))

            cur.execute("DELETE FROM markers " +
                        "WHERE id in %s", (tuple(marker_ids),))

            try:
                for path_to_delete in paths_to_delete:
                    os.unlink(path_to_delete)
            except Exception as e:
                self.ReportWarning(
                    "Exception was raised when deleting markers. Some orphans may remain. Exception message: %s" % e.message)

            self.conn.commit()
            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return False
        except:
            self.conn.rollback()
            raise

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
            table = prettytable.PrettyTable(("marker_id", "name", "description", "size (nt)"))
            table.align = 'l'
            table.hrules = prettytable.FRAME
            table.vrules = prettytable.NONE

            for (_marker_id, name, description, external_id, size) in cur:
                tup = map(str, (external_id, name, description, size))
                table.add_row(tup)
            print table.get_string()

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
                    "WHERE id = %s ", (marker_set_id,))

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
                    "WHERE set_id = %s ", (marker_set_id,))

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

        try:
            gf = GenomeFilter()
            gf.FilterTreeData(self, marker_ids, genome_ids,
                                      comp_threshold, cont_threshold,
                                      taxa_filter,
                                      guaranteed_genome_list_ids, guaranteed_genome_ids,
                                      alignment, individual,
                                      directory, prefix)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        if build_tree:
            pass
            # To Do: should infer tree here

        self.logger.info('Done.')

        return True

    def CalculateMarkersOnProdigalDirAsync(self, marker_ids, prodigal_dir):
        try:

            cur = self.conn.cursor()

            cur.execute("SELECT id, marker_file_location " +
                        "FROM markers " +
                        "WHERE id in %s", (tuple(marker_ids),))

            async_results = {}

            for (marker_id, marker_path) in cur:
                async_results[marker_id] = self.pool.apply_async(
                    MarkerCalculation.CalculateBestMarkerOnProdigalDir,
                    [str(marker_id), marker_path, prodigal_dir]
                )

            return async_results

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def RunProdigalOnGenomeIdAsync(self, genome_id):
        try:
            cur = self.conn.cursor()

            cur.execute("SELECT fasta_file_location " +
                        "FROM genomes " +
                        "WHERE id = %s", (genome_id,))

            (fasta_path,) = cur.fetchone()

            return self.pool.apply_async(
                MarkerCalculation.RunProdigalOnGenomeFasta,
                [fasta_path]
            )

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

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
                    line = line.rstrip()
                    external_ids.append(line)

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

            genome_list_id = self.CreateGenomeListWorking(
                cur, genome_id_list, name, description, owner_id, private)
            if genome_list_id is False:
                raise GenomeDatabaseError("Unable to create new genome list.")

            self.conn.commit()

            return genome_list_id

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # Function: CreateGenomeListWorking
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
    def CreateGenomeListWorking(self, cur, genome_id_list, name, description, owner_id=None, private=None):
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
            cur.execute(
                query, (name, description, owner_id is None, owner_id, private))
            (genome_list_id,) = cur.fetchone()

            query = "INSERT INTO genome_list_contents (list_id, genome_id) VALUES (%s, %s)"
            cur.executemany(query, [(genome_list_id, x)
                                    for x in genome_id_list])

            return genome_list_id
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # Function: GetGenomeIdListFromGenomeListId
    # Given a genome list id, return all the ids of the genomes contained within that genome list.
    #
    # Parameters:
    # genome_list_id - The genome list id of the genome list whose contents needs to be retrieved.
    #
    # Returns:
    # A list of all the genome ids contained within the specified genome list,
    # False on failure.
    def GetGenomeIdListFromGenomeListId(self, genome_list_id):
        try:
            cur = self.conn.cursor()

            cur.execute("SELECT id, owner_id, owned_by_root, private " +
                        "FROM genome_lists " +
                        "WHERE id = %s ", (genome_list_id,))
            result = cur.fetchone()

            if not result:
                raise GenomeDatabaseError(
                    "No genome list with id: %s" % str(genome_list_id))
            else:
                (_list_id, owner_id, owned_by_root, private) = result
                if private and (not self.currentUser.isRootUser()) and (owned_by_root or owner_id != self.currentUser.getUserId()):
                    raise GenomeDatabaseError(
                        "Insufficient permission to view genome list: %s" % str(genome_list_id))

            cur.execute("SELECT genome_id " +
                        "FROM genome_list_contents " +
                        "WHERE list_id = %s", (genome_list_id,))

            return [genome_id for (genome_id,) in cur.fetchall()]

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # Function: GetGenomeIdListFromGenomeListIds
    # Given a list of ids, return all the ids of the genomes contained
    # within that genome list.
    #
    # Parameters:
    #     genome_list_ids - A list of genome list ids whose contents needs to be retrieved.
    #
    # Returns:
    # A list of all the genome ids contained within the specified genome
    # list(s), False on failure.
    def GetGenomeIdListFromGenomeListIds(self, genome_list_ids):
        try:
            cur = self.conn.cursor()

            temp_table_name = Tools.generateTempTableName()

            if genome_list_ids:
                cur.execute("CREATE TEMP TABLE %s (id integer)" %
                            (temp_table_name,))
                query = "INSERT INTO {0} (id) VALUES (%s)".format(
                    temp_table_name)
                cur.executemany(query, [(x,) for x in genome_list_ids])
            else:
                raise GenomeDatabaseError(
                    "No genome lists given. Can not retrieve IDs")

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
                raise GenomeDatabaseError(
                    "Unknown genome list id(s) given. %s" % str(missing_list_ids))

            # Find any genome list ids that we dont have permission to view
            cur.execute("SELECT id, owner_id, owned_by_root, private " +
                        "FROM genome_lists " +
                        "WHERE id in %s ", (tuple(genome_list_ids),))

            no_permission_list_ids = []
            for (list_id, owner_id, owned_by_root, private) in cur:
                if private and (not self.currentUser.isRootUser()) and (owned_by_root or owner_id != self.currentUser.getUserId()):
                    no_permission_list_ids.append(list_id)

            if no_permission_list_ids:
                raise GenomeDatabaseError(
                    "Insufficient permission to view genome lists: %s" % str(no_permission_list_ids))

            cur.execute("SELECT genome_id " +
                        "FROM genome_list_contents " +
                        "WHERE list_id in %s", (tuple(genome_list_ids),))

            return [genome_id for (genome_id,) in cur.fetchall()]

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

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
            genome_id_list = self.GetGenomeIdListFromGenomeListIds(list_ids)

            if genome_id_list is False:
                raise GenomeDatabaseError(
                    "Unable to view genome list. Can not retrieve genomes IDs for lists: %s" % str(list_ids))

            if not self.PrintGenomesDetails(genome_id_list):
                raise GenomeDatabaseError(
                    "Unable to view genome list. Printing to screen failed of genome ids. %s" % str(genome_id_list))

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def PrintGenomeListsDetails(self, genome_list_ids):
        try:
            cur = self.conn.cursor()

            if not genome_list_ids:
                raise GenomeDatabaseError("Unable to print genome details: No genomes given.")

            cur.execute(
                "SELECT lists.id, lists.name, lists.description, lists.owned_by_root, users.username, count(contents.list_id) " +
                "FROM genome_lists as lists " +
                "LEFT OUTER JOIN users ON lists.owner_id = users.id " +
                "JOIN genome_list_contents as contents ON contents.list_id = lists.id " +
                "WHERE lists.id in %s " +
                "GROUP by lists.id, users.username " +
                "ORDER by lists.id asc ", (tuple(genome_list_ids),)
            )

            # print table
            table = prettytable.PrettyTable(("list_id", "name", "description", "owner", "genome_count"))
            table.align = 'l'
            table.hrules = prettytable.FRAME
            table.vrules = prettytable.NONE

            for (list_id, name, description, owned_by_root, username, genome_count) in cur:
                tup = map(str, (list_id, name, (description if description else ''), ("root" if owned_by_root else username), genome_count))
                table.add_row(tup)
            print table.get_string()

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    # True if has permission, False if not. None on error.
    def HasPermissionToViewGenomeList(self, genome_list_id):
        try:
            cur = self.conn.cursor()

            cur.execute("SELECT owner_id, owned_by_root, private " +
                        "FROM genome_lists " +
                        "WHERE id = %s ", (genome_list_id,))

            result = cur.fetchone()

            if not result:
                raise GenomeDatabaseError(
                    "No genome list with id: %s" % str(genome_list_id))

            (owner_id, owned_by_root) = result

            if not self.currentUser.isRootUser():
                if private and (owned_by_root or owner_id != self.currentUser.getUserId()):
                    return False

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return None

    # True if has permission, False if not. None on error.
    def HasPermissionToEditGenomeList(self, genome_list_id):
        try:
            cur = self.conn.cursor()

            cur.execute("SELECT owner_id, owned_by_root " +
                        "FROM genome_lists " +
                        "WHERE id = %s ", (genome_list_id,))

            result = cur.fetchone()

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
                        "Root user editing of other users lists not yet implmented.")

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            self.conn.rollback()
            return None

    # True on success, false on failure.
    def EditGenomeList(self, genome_list_id, batchfile=None, genomes_external_ids=None, operation=None, name=None, description=None, private=None):

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

            if not self.EditGenomeListWorking(cur, genome_list_id, genome_ids, operation, name, description, private):
                raise GenomeDatabaseError(
                    "Unable to edit genome list: %s" % genome_list_id)

            self.conn.commit()
            return True

        except GenomeDatabaseError as e:
            self.conn.rollback()
            self.ReportError(e.message)
            return False

    # True on success, false on failure/error.
    def EditGenomeListWorking(self, cur, genome_list_id, genome_ids=None, operation=None, name=None, description=None, private=None):
        try:
            edit_permission = self.HasPermissionToEditGenomeList(
                genome_list_id)
            if edit_permission is None:
                raise GenomeDatabaseError(
                    "Unable to retrieve genome list id for editing. Offending list id: %s" % genome_list_id)
            if edit_permission is False:
                raise GenomeDatabaseError(
                    "Insufficient permissions to edit this genome list. Offending list id: %s" % genome_list_id)

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
                cur.execute("UPDATE genome_lists SET " + update_query +
                            " WHERE id = %s", params + [genome_list_id])

            temp_table_name = Tools.generateTempTableName()

            if operation is not None:
                if len(genome_ids) == 0:
                    raise GenomeDatabaseError(
                        "No genome ids given to perform '%s' operation." % operation)

                cur.execute("CREATE TEMP TABLE %s (id integer)" %
                            (temp_table_name,))
                query = "INSERT INTO {0} (id) VALUES (%s)".format(
                    temp_table_name)
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
                    raise GenomeDatabaseError(
                        "Unknown genome list edit operation: %s" % operation)

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

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
            for line in fh:
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
            table = prettytable.PrettyTable(("set_id", "name", "description", "owner", "marker_count"))
            table.align = 'l'
            table.hrules = prettytable.FRAME
            table.vrules = prettytable.NONE

            for (set_id, name, description, owned_by_root, username, marker_count) in cur:
                tup = map(str, (set_id, name, description, ("root" if owned_by_root else username), marker_count))
                table.add_row(tup)
            print table.get_string()

            return True

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def GetVisibleMarkerSetsByOwner(self, owner_id=None, all_non_private=False):
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
            if all_non_private:
                privacy_condition = "(private = False OR private is NULL)"

            conditional_query += "AND (" + \
                privacy_condition + " OR owner_id = %s)"
            params.append(self.currentUser.getUserId())

        cur.execute("SELECT id " +
                    "FROM marker_sets " +
                    "WHERE 1 = 1 " +
                    conditional_query, params)

        return [set_id for (set_id,) in cur]

    # Returns list of marker set id. False on failure/error.
    def GetAllVisibleMarkerSetIds(self, all_non_private=False):
        cur = self.conn.cursor()

        conditional_query = ""
        params = []

        if not self.currentUser.isRootUser():
            privacy_condition = "private = False"
            if all_non_private:
                privacy_condition = "(private = False OR private is NULL)"

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
            marker_ids = self.GetMarkerIdListFromMarkerListIds(marker_set_ids)

            if marker_ids is None:
                raise GenomeDatabaseError(
                    "Unable to view marker set. Can not retrieve marker IDs for sets: %s" % str(marker_set_ids))

            if not self.PrintMarkersDetails(marker_ids):
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
