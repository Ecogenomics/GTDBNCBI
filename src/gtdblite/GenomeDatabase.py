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

import os
import pwd
import logging
from multiprocessing import Pool

import prettytable

from biolib.external.fasttree import FastTree

import Config
from UserManager import UserManager
from GenomeDatabaseConnection import GenomeDatabaseConnection
from GenomeManager import GenomeManager
from GenomeListManager import GenomeListManager
from MarkerManager import MarkerManager
from MarkerSetManager import MarkerSetManager
from MetadataManager import MetadataManager
from TreeManager import TreeManager
from Exceptions import GenomeDatabaseError
from AlignedMarkerManager import AlignedMarkerManager
from GenomeRepresentativeManager import GenomeRepresentativeManager


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

    def Login(self, user, login_as_root):
        """Login to database."""

        if not self.conn.IsPostgresConnectionActive():
            raise GenomeDatabaseError(
                "Unable to establish database connection")

        cur = self.conn.cursor()
        user_mngr = UserManager(cur, None)

        if user:
            if not user_mngr.rootLogin(self.GetLinuxUsername()):
                raise GenomeDatabaseError(
                    "Unable to impersonate user %s." % user)

            if not user_mngr.userLogin(user):
                raise GenomeDatabaseError(
                    "Unable to impersonate user %s." % user)

        elif login_as_root:
            if not user_mngr.rootLogin(self.GetLinuxUsername()):
                raise GenomeDatabaseError("Unable to become root user.")
        else:
            if not user_mngr.userLogin(self.GetLinuxUsername()):
                raise GenomeDatabaseError("Database login failed.")

        self.currentUser = user_mngr.currentUser

    def GetLinuxUsername(self):
        return pwd.getpwuid(os.getuid())[0]

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

    # TODO: This should not be here, technically the backend is agnostic so
    # shouldn't assume command line.
    def Confirm(self, msg):
        raw = raw_input(msg + " (y/N): ")
        if raw.upper() == "Y":
            return True
        return False

    def addUser(self, username, role, login_has_root):
        try:
            cur = self.conn.cursor()
            user_mngr = UserManager(cur, self.currentUser)
            user_mngr.addUser(username, role, login_has_root)
            self.conn.commit()
            return True
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def editUser(self, username, role, login_has_root):
        try:
            cur = self.conn.cursor()
            user_mngr = UserManager(cur, self.currentUser)
            user_mngr.editUser(username, role, login_has_root)
            self.conn.commit()
            return True
        except GenomeDatabaseError as e:
            self.conn.rollback()
            self.ReportError(e.message)
            return False

    def AddGenomes(self, batchfile,
                   checkm_file,
                   study_file,
                   modify_genome_list_id=None,
                   new_genome_list_name=None):
        """Add genomes to database.

        Parameters
        ----------
        batchfile : str
            Name of file describing genomes to add.
        checkm_file : str
            Name of file containing CheckM quality estimates.
        study_file : str
            Name of file describing study from which genomes were recovered.

        Returns
        -------
        bool
            True if genomes added without error.
        """
        try:
            self.logger.info('Adding genomes to database.')

            cur = self.conn.cursor()
            genome_list_mngr = GenomeListManager(cur, self.currentUser)

            if modify_genome_list_id is not None:
                if new_genome_list_name is not None:
                    raise GenomeDatabaseError(
                        "Unable to both modify and create genome lists at the same time.")

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
            genome_ids = genome_mngr.addGenomes(checkm_file, batchfile, study_file)

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

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False
        except:
            self.conn.rollback()
            raise

        self.logger.info('Done.')
        return True

    # True if has permission. False if doesn't. None on error.
    def DeleteGenomes(self, batchfile=None, external_ids=None, reason=None):
        '''
        Delete genomes from database.

        :param batchfile:Name of file describing genomes to delete.
        :param external_ids: List of ids describing genomes to delete
        :param reason: Reason of deletion for the listed genomes

        Returns
        -------
        bool
            True if genomes deleted without error.
        '''
        try:
            cur = self.conn.cursor()

            if external_ids is None:
                external_ids = []

            if batchfile:
                for line in open(batchfile, "rb"):
                    if line[0] != '#':
                        external_ids.append(line.rstrip().split('\t')[0])

            genome_mngr = GenomeManager(cur, self.currentUser)
            genome_ids = genome_mngr.externalGenomeIdsToGenomeIds(external_ids)
            genome_mngr.deleteGenomes(batchfile, genome_ids, reason)

            self.conn.commit()
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

    # True on success. False on failure/error.
    def ViewGenomes(self, batchfile=None, external_ids=None):
        try:
            cur = self.conn.cursor()
            genome_mngr = GenomeManager(cur, self.currentUser)

            genome_ids = []
            if external_ids is None and batchfile is None:
                genome_ids = genome_mngr.allGenomeIds()
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

                genome_ids = genome_mngr.externalGenomeIdsToGenomeIds(
                    external_ids)
                if genome_ids is None:
                    raise GenomeDatabaseError("Can not retrieve genome ids.")

            # print genome details
            header, rows = genome_mngr.printGenomeDetails(genome_ids)
            self.PrintTable(header, rows)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

    def ViewMarkers(self, batchfile=None, external_ids=None):
        try:
            cur = self.conn.cursor()

            marker_ids = []
            if external_ids is None and batchfile is None:
                marker_set_mngr = MarkerSetManager(cur, self.currentUser)
                marker_ids = marker_set_mngr.getAllMarkerIds()
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

                marker_mngr = MarkerManager(cur, self.currentUser)
                marker_ids = marker_mngr.externalMarkerIdsToMarkerIds(
                    external_ids)
                if marker_ids is False:
                    raise GenomeDatabaseError("Can not retrieve marker ids.")

            header, rows = marker_mngr.printMarkerDetails(marker_ids)
            self.PrintTable(header, rows)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

    def MakeTreeData(self, marker_ids, genome_ids,
                     directory, prefix,
                     quality_threshold, comp_threshold, cont_threshold,
                     taxa_filter,
                     excluded_genome_list_ids, excluded_genome_ids,
                     guaranteed_genome_list_ids, guaranteed_genome_ids,
                     alignment, individual,
                     build_tree=True):

        try:
            cur = self.conn.cursor()

            # ensure all genomes have been assigned to a representatives
            genome_rep_mngr = GenomeRepresentativeManager(cur, self.currentUser, self.threads)
            genome_rep_mngr.assignToRepresentative()

            # create tree data
            self.logger.info('Creating tree data for %d genomes using %d marker genes.' %
                                (len(genome_ids), len(marker_ids)))

            tree_mngr = TreeManager(cur, self.currentUser)
            genomes_to_retain, chosen_markers_order, chosen_markers = tree_mngr.filterGenomes(marker_ids, genome_ids,
                                                                                              quality_threshold, comp_threshold, cont_threshold,
                                                                                              taxa_filter,
                                                                                              excluded_genome_list_ids,
                                                                                              excluded_genome_ids,
                                                                                              guaranteed_genome_list_ids,
                                                                                              guaranteed_genome_ids)

            aligned_mngr = AlignedMarkerManager(cur, self.threads)
            aligned_mngr.calculateAlignedMarkerSets(genomes_to_retain, marker_ids)

            msa_file = tree_mngr.writeFiles(marker_ids,
                                             genomes_to_retain,
                                             directory,
                                             prefix,
                                             chosen_markers_order,
                                             chosen_markers,
                                             alignment,
                                             individual)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        if build_tree:
            self.logger.info('Inferring tree for %d genome under WAG and GAMMA models.' % len(genomes_to_retain))

            output_tree = os.path.join(directory, prefix + '_phylogeny.wag_gamma.tree')
            output_tree_log = os.path.join(directory, prefix + '_fasttree.log')
            log_file = os.path.join(directory, prefix + '_fasttree_output.txt')

            fasttree = FastTree(multithreaded=True)
            fasttree.run(msa_file, 'prot', 'wag', output_tree, output_tree_log, log_file)

        self.logger.info('Done.')

        return True

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

            genome_mngr = GenomeManager(cur, self.currentUser)
            genome_id_list = genome_mngr.externalGenomeIdsToGenomeIds(
                external_ids)
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
            genome_id_list = genome_list_mngr.getGenomeIdsFromGenomeListIds(list_ids)

            if not genome_id_list:
                raise GenomeDatabaseError(
                    "Unable to view genome list. Cannot retrieve genomes IDs for lists: %s" % str(list_ids))

            genome_mngr = GenomeManager(cur, self.currentUser)
            header, rows = genome_mngr.printGenomeDetails(genome_id_list)
            self.PrintTable(header, rows)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)

        return True

    def PrintGenomeListsDetails(self, genome_list_ids):
        """Print genome list details.

        Parameters
        ----------
        genome_list_ids : iterable
            Unique identifier of genome lists in database.

        Returns
        -------
        bool
            True if successful, else False.
        """

        try:
            cur = self.conn.cursor()

            genome_list_mngr = GenomeListManager(cur, self.currentUser)
            header, rows = genome_list_mngr.printGenomeListsDetails(
                genome_list_ids)

            self.PrintTable(header, rows)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

    def printMarkerSetsDetails(self, marker_set_ids):
        """Print genome list details.

        Parameters
        ----------
        marker_set_ids : iterable
            Unique identifier of marker sets in database.

        Returns
        -------
        bool
            True if successful, else False.
        """

        try:
            cur = self.conn.cursor()

            marker_set_mngr = MarkerSetManager(cur, self.currentUser)
            header, rows = marker_set_mngr.printMarkerSetsDetails(
                marker_set_ids)

            self.PrintTable(header, rows)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

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
                genome_mngr = GenomeManager(cur, self.currentUser)
                genome_ids = genome_mngr.externalGenomeIdsToGenomeIds(
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

    def DeleteGenomeLists(self, list_ids=None):
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

            marker_mngr = MarkerManager(cur, self.currentUser)
            marker_id_list = marker_mngr.externalMarkerIdsToMarkerIds(
                external_ids)
            if marker_id_list is False:
                raise GenomeDatabaseError(
                    "Unable to retreive marker ids for provided markers.")

            owner_id = None
            if not self.currentUser.isRootUser():
                owner_id = self.currentUser.getUserId()

            marker_set_mngr = MarkerSetManager(cur, self.currentUser)
            marker_list_id = marker_set_mngr.createMarkerSet(
                marker_id_list, name, description, owner_id, private)
            if marker_list_id is False:
                raise GenomeDatabaseError("Unable to create new marker set.")

# TODEL            self.PrintMarkerSetsDetails(marker_list_id)

            self.conn.commit()

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return marker_list_id

    def EditMarkerSet(self, marker_set_id, batchfile=None, marker_external_ids=None, operation=None, name=None, description=None, private=None):
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
            cur = self.conn.cursor()

            if batchfile:
                if marker_external_ids is None:
                    marker_external_ids = []
                for line in batchfile:
                    line = line.rstrip()
                    marker_external_ids.append(line)

            if marker_external_ids is not None:
                marker_set_mngr = MarkerManager(cur, self.currentUser)
                marker_external_ids = marker_set_mngr.externalMarkerIdsToMarkerIds(
                    marker_external_ids)

            marker_set_mngr = MarkerSetManager(cur, self.currentUser)
            if not marker_set_mngr.editMarkerSet(marker_set_id, marker_external_ids, operation, name, description, private):
                raise GenomeDatabaseError(
                    "Unable to edit marker set: %s" % marker_set_id)

            self.conn.commit()
        except GenomeDatabaseError as e:
            self.conn.rollback()
            self.ReportError(e.message)
            return False

        return True

    def deleteMarkerSets(self, set_ids=None):
        '''
        Delete a list of marker sets

        :param list_ids: List of marker IDs to delete
        '''
        try:
            cur = self.conn.cursor()
            marker_set_mngr = MarkerSetManager(cur, self.currentUser)
            if not marker_set_mngr.deleteMarkerSets(set_ids):
                raise GenomeDatabaseError(
                    "Unable to edit marker sets")
            self.conn.commit()
        except GenomeDatabaseError as e:
            self.conn.rollback()
            self.ReportError(e.message)
            return False
        return True

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

    def ViewMarkerSetsContents(self, marker_set_ids):
        try:
            cur = self.conn.cursor()

            marker_set_mngr = MarkerSetManager(cur, self.currentUser)
            marker_ids = marker_set_mngr.getMarkerIdListFromMarkerSetIds(marker_set_ids)

            if marker_ids is None:
                raise GenomeDatabaseError(
                    "Unable to view marker set. Can not retrieve marker IDs for sets: %s" % str(marker_set_ids))

            marker_manager = MarkerManager(cur, self.currentUser)
            header, rows = marker_manager.printMarkerDetails(marker_ids)
            self.PrintTable(header, rows)

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

    def ViewMetadata(self):
        try:
            cur = self.conn.cursor()
            metaman = MetadataManager(cur, self.currentUser)
            metaman.viewMetadata()
            return True
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def ExportMetadata(self, path):
        try:
            cur = self.conn.cursor()

            # ensure all genomes have been assigned to a representatives
            genome_rep_mngr = GenomeRepresentativeManager(cur, self.currentUser, self.threads)
            genome_rep_mngr.assignToRepresentative()

            metaman = MetadataManager(cur, self.currentUser)
            metaman.exportMetadata(path)

            self.conn.commit()
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

    def ImportMetadata(self, table=None, field=None, typemeta=None, metafile=None):
        try:
            cur = self.conn.cursor()
            metaman = MetadataManager(cur, self.currentUser)
            metaman.importMetadata(table, field, typemeta, metafile)
            self.conn.commit()
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True

    def CreateMetadata(self, path):
        try:
            cur = self.conn.cursor()
            metaman = MetadataManager(cur, self.currentUser)
            metaman.createMetadata(path)
            self.conn.commit()
            return True
        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def ReportStats(self):
        """Report general database statistics."""

        try:
            cur = self.conn.cursor()

            # make sure representative have been determine for all genomes
            genome_rep_mngr = GenomeRepresentativeManager(cur, self.currentUser, self.threads)
            genome_rep_mngr.assignToRepresentative()

            # determine number of genomes from different database sources
            cur.execute("SELECT name, external_id_prefix, " +
                        "(SELECT COUNT(*) "
                        "FROM genomes "
                        "WHERE genomes.genome_source_id = genome_sources.id) "
                        "FROM genome_sources " +
                        "ORDER BY id")
            genome_counts = cur.fetchall()

            print ''
            header = ('Genome Source', 'Prefix', 'Genome Count')
            rows = []
            for tup in genome_counts:
                rows.append(tup)

            self.PrintTable(header, genome_counts)

            print ''
            print 'Total genomes: %d' % (sum([x[2] for x in genome_counts]))

            # report number of reference genomes
            cur.execute("SELECT COUNT(*) "
                        "FROM metadata_taxonomy "
                        "WHERE gtdb_representative = True")
            ref_genome_counts = cur.fetchone()[0]
            print 'Number of reference genomes: %d' % ref_genome_counts

            self.conn.commit()

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

        return True
