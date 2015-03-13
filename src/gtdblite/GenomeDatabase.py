import hashlib
import shutil
import os
import datetime
import time

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
            self.ReportError("User not found")
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
        
        cur.execute("SELECT genome_source_id, id_at_source, tree_id_prefix " +
                    "FROM genomes, genome_sources " +
                    "WHERE id = %s "+
                    "AND genome_source_id = genome_sources.id", (genome_id,))
        
        result = cur.fetchall()
        if len(result) == 0:
            self.ReportError("Genome id not found: %s." % genome_id)
            return False
        
        for (genome_source_id, id_at_source, tree_id_prefix) in result:
            target_file_name = tree_id_prefix + "_" + str(id_at_source)
            try:
                shutil.copy(fasta_file, os.path.join(self.genomeCopyDir, target_file_name))
            except:
                self.ReportError("Copy to genome storage dir failed.")
                return False
            return target_file_name
        
      
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
            
            if not self.ModifyGenomeListWorking(cur, modify_genome_list_id, genome_ids=added_genome_ids, operation='add'):
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
                fasta_fh.close()
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
                
            cur.execute("SELECT id, tree_id_prefix, user_accessible FROM genome_sources WHERE name = %s" , (source,))
            source_id = None
            prefix = None
            
            for (id, tree_id_prefix, user_accessible) in cur:
                if (not user_accessible):
                    if id_at_source == None:
                        raise GenomeDatabaseError("Cannot auto generate ids at source for the %s genome source." % source)
                    if (not self.currentUser.isRootUser()):
                        raise GenomeDatabaseError("Only the root user can add genomes to the %s genome source." % source)
                source_id = id
                prefix = tree_id_prefix
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

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return None
        except:            
            raise 
        
        genome_id = cur.fetchone()[0]
        return genome_id
        
    
    def AddMarkers(self, batchfile, force_overwrite=False):
    
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

    def AddMarkersWorking(self, cur, marker_file_path, name, desc, force_overwrite=False,
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
                
        cur.execute("SELECT marker_id " +
                    "FROM marker_set_contents " +
                    "WHERE set_id = %s ", (marker_set_id,))
        
        marker_id_list = [x[0] for x in cur.fetchall()]
        
        return marker_id_list
  
        
    def MakeTreeData(self, marker_list, list_of_genome_ids, directory, prefix, profile=None, config_dict=None, build_tree=True):    
       
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
        
        for genome_id in list_of_genome_ids:
            uncalculated = self.FindUncalculatedMarkersForGenomeId(genome_id, marker_list)
            if len(uncalculated) != 0:
                uncalculated_marker_dict[genome_id] = uncalculated
                uncalculated_marker_count += len(uncalculated)
        
        print "%i genomes contain %i uncalculated markers." % (len(uncalculated_marker_dict.keys()), uncalculated_marker_count)
        # TODO: Add Confirm step

        for (genome_id, uncalculated) in incalculated_marker_dict:
            self.RecalculateMarkersForGenome(genome_id, uncalculated)
        
        return profiles.profiles[profile].MakeTreeData(self, marker_set_id, list_of_genome_ids,
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
    #     genome_list_id - The genome list id of the genome list whose contents needs to be retrieved.
    #
    # Returns:
    #   A list of all the genome ids contained within the specified genome list, None on failure.
    def GetGenomeIdListFromGenomeListId(self, genome_list_id):
        
        cur = self.conn.cursor()
        
        cur.execute("SELECT id " +
                    "FROM genome_lists " +
                    "WHERE id = %s", (genome_list_id,))
        
        if not cur.fetchone():
            self.ReportError("No genome list with id: " + str(genome_list_id))
            return None
        
        cur.execute("SELECT genome_id " +
                    "FROM genome_list_contents " +
                    "WHERE list_id = %s", (genome_list_id,))
        
        result = cur.fetchall()
        
        return [genome_id for (genome_id,) in result]
    
    
    # Function: ModifyGenomeList
    # Modify the details or contents of an existing genome list.
    #
    # Parameters:
    #     genome_list_id - The genome list id of the genome list to modify.
    #     name - If not None, update the genome list's name to this.
    #     description - If not None, update the genome list's description to this.
    #     genome_ids - List of genome id that will modify the contents of the genome list (see operation parameter)
    #     operation - Perform this operation on the genome ids given in the genome_ids parameter with respect to the genome list (options: add, remove)
    #     private - If not True, change if this list is private.
    #
    # Returns:
    #   True on success, False otherwise.
    def ModifyGenomeListWorking(self, cur, genome_list_id, name=None, description=None, genome_ids=None,
                                operation=None, private=None):
        
        query = "SELECT owned_by_root, owner_id FROM genome_lists WHERE id = %s";
        cur.execute(query, (genome_list_id,))
        result = cur.fetchone()
        if not result:
            self.ReportError("Can't find specified genome list Id: " + str(genome_list_id))
            return False
        
        (owned_by_root, owner_id) = result
        
        # Need to check permissions to edit this list.
        if self.currentUser.isRootUser():
            if not owned_by_root:
                self.ReportError("Root user editing of other users lists not yet implmented.")
                return False
        else:
            if owned_by_root:
                self.ReportError("Only the root user can edit root owned lists.")
                return False
            if owner_id != self.currentUser.getUserId():
                self.ReportError("Insufficient privileges to edit this genome list.")
                return False
        
        if name is not None:
            query = "UPDATE genome_lists SET name = %s WHERE id = %s";
            cur.execute(query, (name, genome_list_id))
            
        if description is not None:
            query = "UPDATE genome_lists SET description = %s WHERE id = %s";
            cur.execute(query, (description, genome_list_id))
            
        if private is not None:
            query = "UPDATE genome_lists SET private = %s WHERE id = %s";
            cur.execute(query, (private, genome_list_id))
            
        temp_table_name = "TEMP" + str(int(time.time()))
        
        if genome_ids:
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
        
        return True
    
    
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
            conditional_query += "AND (list.private = False OR owner_id = %s)"
            params.append(self.currentUser.getUserId())
        
        cur.execute(
            "SELECT list.id, list.name, list.description, users.username, count(contents.list_id) " +
            "FROM genome_lists as list " +
            "LEFT OUTER JOIN users ON list.owner_id = users.id " + 
            "JOIN genome_list_contents as contents ON contents.list_id = list.id " +
            "WHERE 1 = 1 " +
            conditional_query +
            "GROUP by list.id, users.username " +
            "ORDER by list.id asc " ,
            params
        )
        
        print cur.query
        
        return cur.fetchall()
    
    
    def GetAllVisibleGenomeLists(self):
        cur = self.conn.cursor()

        conditional_query = ""
        params = []
            
        if not self.currentUser.isRootUser():
            conditional_query += "AND (list.private = False OR owner_id = %s)"
            params.append(self.currentUser.getUserId())
        
        cur.execute(
            "SELECT list.id, list.name, list.description, users.username, count(contents.list_id) " +
            "FROM genome_lists as list " +
            "LEFT OUTER JOIN users ON list.owner_id = users.id " + 
            "JOIN genome_list_contents as contents ON contents.list_id = list.id " +
            "WHERE 1 = 1 " +
            conditional_query +
            "GROUP by list.id, users.username " +
            "ORDER by list.id asc " ,
            params
        )
        
        print cur.query
        
        return cur.fetchall()
    
    
    
    
    
    
    
    
    
    