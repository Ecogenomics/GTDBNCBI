import hashlib
import shutil
import os
import datetime

from gtdblite.User import User
from gtdblite.GenomeDatabaseConnection import GenomeDatabaseConnection

class GenomeDatabase(object):
    def __init__(self):
        self.conn = GenomeDatabaseConnection()
        self.currentUser = None
        self.lastErrorMessage = None
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
        self.lastErrorMessage = str(msg) + "\n"
    
        
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
        
      
    def AddFastaGenome(self, fasta_file, copy_fasta, name, desc, force_overwrite=False, source=None, id_at_source=None, completeness=0, contamination=0):
        cur = self.conn.cursor()
        
        genome_id = self.AddFastaGenomeWorking(cur, fasta_file, name, desc, force_overwrite, source, id_at_source, completeness, contamination)
        if (genome_id):
            if copy_fasta:
                target_file_name = self.CopyFastaToCopyDir(fasta_file, genome_id)
                if target_file_name:
                    self.conn.commit()
                    return True
            else:
                self.conn.commit()
                return True
        
        self.conn.rollback()
        return False        
    
    def AddManyFastaGenomes(self, batchfile, checkM_file, force_overwrite=False):
    
        checkm_fh = open(checkM_file, "rb")
        
        expected_headers = ["Bin Id", "Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2",
                            "3", "4", "5+", "Completeness", "Contamination",  "Strain heterogeneity"]
        
        # Check the CheckM headers are consistent
        if checkm_fh.readline().rstrip() != "\t".join(expected_headers):
            self.ReportError("CheckM file header inconsistent. Expected: \n\t%s\nGot:\n\t%s\n" %
                             ("\t".join(expected_headers), checkm_fh.readline().rstrip()))
            return False
        
        # Populate CheckM results dict
        checkM_results_dict = {}
        
        for line in checkm_fh:
            line = line.rstrip()
            splitline = line.split("\t")
            file_name, completeness, contamination = splitline[0], splitline[11], splitline[12]

            checkM_results_dict[file_name] = {"completeness" : completeness, "contamination" : contamination}
        
        checkm_fh.close()
        
        # Add the genomes
        cur = self.conn.cursor()
        fasta_paths_to_copy = []
        
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
                self.ReportError("Couldn't find checkM result for %s (%s)" % (name,abs_path))
                return False
            
            genome_id = self.AddFastaGenomeWorking(
                cur, abs_path, name, desc, force_overwrite, source_name, id_at_source,
                checkM_results_dict[basename]["completeness"], checkM_results_dict[basename]["contamination"]
            )
            
            # Rollback everything if addition fails
            if not (genome_id):
                self.conn.rollback()
                return False
            
            fasta_paths_to_copy.append(abs_path)
    
        #if copy_fastas:
        #    # TODO: Copy the fastas if required, rollback if fails
        #    pass
        
        self.conn.commit()
        return True
    
    def AddFastaGenomeWorking(self, cur, fasta_file_path, name, desc, force_overwrite=False,
                              source=None, id_at_source=None, completeness=0, contamination=0):
        try:
            fasta_fh = open(fasta_file_path, "rb")
        except:
            self.ReportError("Cannot open Fasta file: " + fasta_file_path)
            fasta_fh.close()
            return None
        
        m = hashlib.sha256()
        for line in fasta_fh:
            m.update(line)                
        fasta_sha256_checksum = m.hexdigest()
        fasta_fh.close()

        if source is None:
            source = self.defaultGenomeSourceName
        
        cur.execute("SELECT id, tree_id_prefix, user_accessible FROM genome_sources WHERE name = %s" , (source,))
        source_id = None
        prefix = None
        
        for (id, tree_id_prefix, user_accessible) in cur:
            if (not user_accessible):
                if id_at_source == None:
                    self.ReportError("Cannot auto generate ids at source for the %s genome source." % source)
                    return None
                if (not self.currentUser.isRootUser()):
                    self.ReportError("Only the root user can add genomes to the %s genome source." % source)
                    return None
            source_id = id
            prefix = tree_id_prefix
            break
        
        if source_id is None:
            self.ReportError("Could not find the %s genome source." % source)
            return None
        
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
                self.ReportError("Force overwrite not implemented yet")
                return None
            else:
                self.ReportError("Genome source '%s' already contains id '%s'. Use -f to force an overwrite." % (source, id_at_source))
                return None
        else:
            cur.execute("INSERT INTO genomes " + columns + " "
                        "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) " + 
                        "RETURNING id" ,
                        (name, desc, self.currentUser.isRootUser(), owner_id, fasta_file_path, fasta_sha256_checksum, source_id, id_at_source, added, completeness, contamination))
        
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
        
    def MakeTreeData(self, marker_list, list_of_genome_ids, directory, profile=None, prefix=None, config_dict=None):    
       
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    