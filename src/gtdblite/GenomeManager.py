'''
Created on 21 Oct 2015

@author: la.p.chaumeil

'''

import os
import logging


from gtdblite.GenomeDatabaseConnection import GenomeDatabaseConnection
from gtdblite.Exceptions import GenomeDatabaseError
from Tools import fastaPathGenerator, generateTempTableName
from gtdblite import Config

from biolib.common import make_sure_path_exists


class GenomeManager(object):
    '''
    classdocs
    '''

    def __init__(self, user):
        '''
        Constructor
        '''
        self.conn = GenomeDatabaseConnection()
        self.conn.MakePostgresConnection()
        self.currentUser = user
        self.logger = logging.getLogger()

        self.deprecatedUserDir = None
        if Config.GTDB_DPRCTD_USR_DIR:
            self.deprecatedUserDir = Config.GTDB_DPRCTD_USR_DIR

        self.deprecatedNCBIDir = None
        if Config.GTDB_DPRCTD_NCBI_DIR:
            self.deprecatedNCBIDir = Config.GTDB_DPRCTD_NCBI_DIR

    # TODO: This should not be here, techincally the backend is agnostic so
    # shouldn't assume command line.
    def Confirm(self, msg):
        raw = raw_input(msg + " (y/N): ")
        if raw.upper() == "Y":
            return True
        return False

    def loggerSetup(self, silent=False):
        """Set logging for application.

        Parameters
        ----------
        output_dir : str
            Output directory for log file.
        silent : boolean
            Flag indicating if output to stdout should be suppressed.
        """

        # setup general properties of logger
        logger = logging.getLogger('')
        logger.setLevel(logging.DEBUG)
        log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                       datefmt="%Y-%m-%d %H:%M:%S")

        file_logger = logging.FileHandler(os.path.join(
            os.path.realpath(__file__), '..', '..', '..', 'logs', 'genome_manager.log'), 'a')
        file_logger.setFormatter(log_format)
        logger.addHandler(file_logger)

    def DeleteGenomes(self, batchfile=None, external_ids=None, reason=None):
        self.loggerSetup()
        '''
        Delete Genomes
        Returns True for success or False for fail

        Parameters:
        :param batchfile: text file listing a range of ids to delete
        :param external_ids: a list of ids can be written directly in the command line
        '''
        try:
            cur = self.conn.cursor()

            if external_ids is None:
                external_ids = []

            if batchfile:
                fh = open(batchfile, "rb")
                external_ids.extend([line.rstrip('\n') for line in fh])
                fh.close()

            genome_ids = self.ExternalGenomeIdsToGenomeIds(external_ids)

            if genome_ids is False:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Unable to retrieve genome ids.")

            has_permission, username, genomes_owners = self.HasPermissionToEditGenomes(
                genome_ids)

            print genome_ids

            if has_permission is None:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Unable to retrieve permissions for genomes.")

            if has_permission is False:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Insufficient permissions.")

            if not self.Confirm("Are you sure you want to delete %i genomes (this action cannot be undone)" % len(genome_ids)):
                raise GenomeDatabaseError("User aborted database action.")

            cur.execute("DELETE FROM aligned_markers " +
                        "WHERE genome_id in %s ", (tuple(genome_ids),))

            cur.execute("DELETE FROM genome_list_contents " +
                        "WHERE genome_id in %s", (tuple(genome_ids),))

            # Deletion of metadata

            cur.execute("DELETE FROM metadata_genes " +
                        "WHERE id in %s", (tuple(genome_ids),))
            cur.execute("DELETE FROM metadata_ncbi " +
                        "WHERE id in %s", (tuple(genome_ids),))
            cur.execute("DELETE FROM metadata_nucleotide " +
                        "WHERE id in %s", (tuple(genome_ids),))
            cur.execute("DELETE FROM metadata_taxonomy " +
                        "WHERE id in %s", (tuple(genome_ids),))

            cur.execute("DELETE FROM genomes " +
                        "WHERE id in %s", (tuple(genome_ids),))

            for genome, info in genomes_owners.iteritems():
                if str(username) != str(info.get("owner")):
                    logging.info('''Genome {0} has been deleted by {1} for the following reason '{2}'
                                      WARNING: {1} is not the owner of this {0} (real owner {3} )
                                      {0} needs to be moved manually to the deprecated folder'''.format(genome, username, reason, info.get("owner")))
                else:
                    if info.get("prefix") is "U":
                        target = os.path.dirname(
                            os.path.join(self.deprecatedUserDir, info.get("relative_path")))
                    elif info.get("prefix") is "NCBI":
                        target = os.path.join(
                            self.deprecatedNCBIDir, info.get("relative_path"))
                    make_sure_path_exists(target)
                    os.rename(
                        os.path.dirname(fastaPathGenerator(info.get("relative_path"), info.get("prefix"))), target)
                    logging.info("Genome {0} has been deleted by {1} for the following reason '{2}'".format(
                        genome, username, reason))

            self.conn.commit()
            return True

        except GenomeDatabaseError as e:
            self.conn.rollback()
            raise e

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

            temp_table_name = generateTempTableName()

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
                temp_table_name = generateTempTableName()
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

    # True if has permission. False if doesn't. None on error.
    def HasPermissionToEditGenomes(self, genome_ids):
        '''
        Function HasPermissionToEditGenomes
        Check if a user is entitled to delete genomes.
        Users can delete their own genomes, Admin can delete any genomes

        :param genome_ids:list of genomes is to delete

        Return a tuple containing:
        - Boolean returning the state of the function
        - The username currently running the delete function
        -a Dictionary listing the list of genomes and where each genome has saved 
            the owner,the Prefix (U or NCBI) , the relative path 

        '''
        try:
            cur = self.conn.cursor()

            if not genome_ids:
                raise GenomeDatabaseError(
                    "Unable to retrieve genome permissions, no genomes given: %s" % str(genome_ids))

            cur.execute("SELECT gs.external_id_prefix,gs.external_id_prefix || '_'|| genomes.id_at_source, owner_id, username, owned_by_root,fasta_file_location "
                        "FROM genomes " +
                        "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
                        "LEFT JOIN genome_sources gs ON gs.id = genomes.genome_source_id " +
                        "WHERE genomes.id in %s", (tuple(genome_ids),))

            dict_genomes_user = {}
            for (prefix, public_id, owner_id, username, owned_by_root, fasta_path) in cur:

                if not self.currentUser.isRootUser():
                    if (owned_by_root or owner_id != self.currentUser.getUserId()):
                        print (
                            "WARNING: Insufficient permissions to edit genome {0}".format(public_id))
                        print logging.warn("{0} is trying to delete genome {1} owned by {2}".format(self.currentUser.getUsername(), public_id, username))
                        return (False, None, None)
                dict_genomes_user[public_id] = {
                    "owner": username, "prefix": prefix, "relative_path": fasta_path}

            if self.currentUser.isRootUser():
                current_username = self.currentUser.getElevatedFromUsername()
            else:
                current_username = self.currentUser.getUsername()

            return (True, current_username, dict_genomes_user)

        except GenomeDatabaseError as e:
            raise e
