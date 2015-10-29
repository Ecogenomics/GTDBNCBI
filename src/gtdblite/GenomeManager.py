<<<<<<< HEAD
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
import tempfile
import shutil
import logging
import datetime

import psycopg2
from psycopg2.extensions import AsIs

from gtdblite import Config
from gtdblite import ConfigMetadata
from gtdblite.GenomeDatabaseConnection import GenomeDatabaseConnection
from gtdblite.Exceptions import GenomeDatabaseError
from gtdblite.MetadataManager import MetadataManager
from gtdblite import Tools
from gtdblite.Prodigal import Prodigal
from gtdblite.TigrfamSearch import TigrfamSearch
from gtdblite.PfamSearch import PfamSearch

from biolib.checksum import sha256
from Tools import fastaPathGenerator, generateTempTableName
from biolib.common import make_sure_path_exists


class GenomeManager(object):
    """Manage the processing of new genomes and querying genome information."""

    def __init__(self, currentUser, threads):
        """Initialize.

        Parameters
        ----------
        currentUser : User
            Current user of database.
        threads : int
            Number of threads to use for processing.
        """

        self.logger = logging.getLogger()

        self.currentUser = currentUser
        self.threads = threads

        self.conn = GenomeDatabaseConnection()
        self.conn.MakePostgresConnection()

        self.defaultGenomeSourceName = 'user'

        self.genomeFileSuffix = ConfigMetadata.GENOME_FILE_SUFFIX
        self.proteinFileSuffix = ConfigMetadata.PROTEIN_FILE_SUFFIX
        self.ntGeneFileSuffix = ConfigMetadata.NT_GENE_FILE_SUFFIX
        self.gffFileSuffix = ConfigMetadata.GFF_FILE_SUFFIX
        self.tranTableFileSuffix = ConfigMetadata.TRANSLATION_TABLE_SUFFIX
        self.checksumSuffix = ConfigMetadata.CHECKSUM_SUFFIX

        self.genomeCopyDir = Config.GTDB_GENOME_COPY_DIR
        self.deprecatedUserDir = Config.GTDB_DPRCTD_USR_DIR
        self.deprecatedNCBIDir = Config.GTDB_DPRCTD_NCBI_DIR
        
    def _loggerSetup(self, silent=False):
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


    def addGenomes(self, cur, checkm_file, batchfile):
        """Add new genomes to DB.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        checkm_file : str
            Name of file containing CheckM results.
        batchfile : str
            Name of file describing genomes to add.

        Returns
        -------
        list
            List of database genome identifiers of added genomes.
        """

        try:
            tmp_output_dir = tempfile.mkdtemp()

            self.logger.info("Reading CheckM file.")
            checkm_results_dict = self._processCheckM(checkm_file)

            genomic_files = self._addGenomeBatch(cur, batchfile, tmp_output_dir)

            self.logger.info("Running Prodigal to identify genes.")
            prodigal = Prodigal(self.threads)
            file_paths = prodigal.run(genomic_files)

            self.logger.info("Calculating and storing metadata for each genomes.")
            metadata_mngr = MetadataManager()
            for db_genome_id, values in genomic_files.iteritems():
                genome_file_paths = file_paths[db_genome_id]
                output_dir, _file = os.path.split(genome_file_paths["aa_gene_path"])

                bin_id = values['checkm_bin_id']
                if bin_id not in checkm_results_dict:
                    raise GenomeDatabaseError("Couldn't find CheckM result for bin  %s." % bin_id)

                metadata_mngr.addMetadata(cur,
                                          db_genome_id,
                                          genome_file_paths["fasta_path"],
                                          genome_file_paths["gff_path"],
                                          checkm_results_dict[bin_id],
                                          output_dir)

            self.logger.info("Identifying TIGRfam protein families.")
            gene_files = [file_paths[db_genome_id]['aa_gene_path'] for db_genome_id in genomic_files]
            tigr_search = TigrfamSearch(self.threads)
            tigr_search.run(gene_files)

            self.logger.info("Identifying Pfam protein families.")
            pfam_search = PfamSearch(self.threads)
            pfam_search.run(gene_files)

            # all genomes were process successfully so move them
            # into the GTDB directory structure
            self.logger.info("Moving files to GTDB directory structure.")
            self._moveGenomes(cur, genomic_files.keys(), tmp_output_dir)
        except:
            if os.path.exists(tmp_output_dir):
                shutil.rmtree(tmp_output_dir)
            raise

        return genomic_files.keys()

    def _moveGenomes(self, cur, db_genome_ids, tmp_output_dir):
        """Move genome files into database directory structure.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        db_genome_ids : list
            Unique database identifiers for genomes.
        tmp_output_dir : str
            Temporary directory with genome data to move into database directory structure.
        """

        # get database genome identifiers
        cur.execute("SELECT genomes.id,user_editable, external_id_prefix || '_' || id_at_source as external_id " +
                    "FROM genomes, genome_sources " +
                    "WHERE genome_source_id = genome_sources.id " +
                    "AND genomes.id in %s", (tuple(db_genome_ids),))

        external_id_dict = {}
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
                raise GenomeDatabaseError("Unable to determine user to add genomes under.")

        gtdb_target_dir = os.path.join(self.genomeCopyDir, username)
        for db_genome_id, external_id in external_id_dict.items():
            tmp_genome_dir = os.path.join(tmp_output_dir, external_id)

            genome_target_dir = os.path.join(gtdb_target_dir, external_id)
            if os.path.exists(genome_target_dir):
                raise GenomeDatabaseError("Genome directory already exists: %s" % genome_target_dir)

            shutil.move(tmp_genome_dir, gtdb_target_dir)

            cur.execute("UPDATE genomes SET fasta_file_location = %s , genes_file_location = %s WHERE id = %s", (
                    os.path.join(genome_target_dir, external_id + self.genomeFileSuffix),
                    os.path.join(genome_target_dir, external_id + self.proteinFileSuffix),
                    db_genome_id))

        shutil.rmtree(tmp_output_dir)

    def _addGenomeBatch(self, cur, batchfile, output_dir):
        """Add genomes specific in batch file to DB.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        batchfile : str
            Name of file describing genomes to add.
        output_dir : str
            Output directory.

        Returns
        -------
        dict
            Dictionary indicating the genomic and gene file for each genome.
        """

        # Add the genomes
        genomic_files = {}
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

            (fasta_path, name, desc, gene_path, source_name, id_at_source) = splitline

            if fasta_path is None or fasta_path == '':
                raise GenomeDatabaseError(
                    "Each line in the batch file must specify a path to the genome's fasta file.")

            if name is None or name == '':
                raise GenomeDatabaseError(
                    "Each line in the batch file must specify a name for the genome.")

            abs_fasta_path = os.path.abspath(fasta_path)

            abs_gene_path = None
            if gene_path is not None and gene_path != '':
                abs_gene_path = os.path.abspath(gene_path)

            genome_id = self._addGenomeToDB(cur, abs_fasta_path, name, desc, source_name, id_at_source, abs_gene_path)
            if not (genome_id):
                raise GenomeDatabaseError("Failed to add genome: %s" % abs_fasta_path)

            cur.execute("SELECT external_id_prefix || '_' || id_at_source as external_id " +
                "FROM genomes, genome_sources " +
                "WHERE genome_source_id = genome_sources.id " +
                "AND genomes.id = %s", (genome_id,))

            external_genome_id = cur.fetchone()[0]

            genome_output_dir = os.path.join(output_dir, external_genome_id)
            if not os.path.exists(genome_output_dir):
                os.makedirs(genome_output_dir)

            fasta_target_file = os.path.join(genome_output_dir, external_genome_id + self.genomeFileSuffix)
            shutil.copy(abs_fasta_path, fasta_target_file)

            genes_target_file = None
            if abs_gene_path:
                genes_target_file = os.path.join(output_dir, external_genome_id, external_genome_id + self.proteinFileSuffix)
                shutil.copy(abs_gene_path, genes_target_file)

            genomic_files[genome_id] = {"checkm_bin_id": os.path.splitext(os.path.basename(abs_fasta_path))[0],
                                            "aa_gene_path": genes_target_file,
                                            "fasta_path": fasta_target_file}
        return genomic_files

    def _addGenomeToDB(self, cur, fasta_file_path, name, desc,
                              source, id_at_source, gene_path):
        """Add genome to database.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        fasta_file_path : str
            Path to genome FASTA file with nucleotide sequences.
        name : str
            Desired name of genome.
        desc : str
            Description of genome.
        source : str
            Source of genome.
        id_at_source : int
            ?
        gene_path : str
            Path to called genes in amino acid space.

        Returns
        -------
        str
            Database identifier of genome.
        """
        try:
            fasta_sha256_checksum = sha256(fasta_file_path)

            gene_sha256_checksum = None
            if gene_path is not None:
                gene_sha256_checksum = sha256(gene_path)
            if source is None:
                source = self.defaultGenomeSourceName

            cur.execute("SELECT id, external_id_prefix, user_editable FROM genome_sources WHERE name = %s", (source,))
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
                raise GenomeDatabaseError(
                    "Genome source '%s' already contains id '%s'. Use -f to force an overwrite." % (source, id_at_source))

            cur.execute("INSERT INTO genomes " + columns + " "
                        "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) " +
                        "RETURNING id",
                        (name, desc, self.currentUser.isRootUser(), owner_id, fasta_file_path, fasta_sha256_checksum, gene_path, gene_sha256_checksum, source_id, id_at_source, added))
            (db_genome_id,) = cur.fetchone()

            return db_genome_id

        except GenomeDatabaseError as e:
            self.ReportError(e.message)
            return False

    def _identifyHeadersCheckM(self, checkm_fh):
        """Parse header information from CheckM file.

        Parameters
        ----------
        checkm_fh : file
            Handle to CheckM file.
        """

        required_headers = {
            "Bin Id": None,
            "Completeness": None,
            "Contamination": None,
            "Marker lineage": None,
            "# genomes": None,
            "# markers": None,
            "# marker sets": None,
            "Strain heterogeneity": None
        }

        # Check the CheckM headers are consistent
        split_headers = checkm_fh.readline().rstrip().split("\t")
        for pos in range(0, len(split_headers)):
            header = split_headers[pos]
            if header not in required_headers:
                continue

            if required_headers[header] is not None:
                raise GenomeDatabaseError(
                    "Seen %s header twice in the CheckM file. Check that the CheckM file is correct: %s." % (header, checkm_fh.name))

            required_headers[header] = pos

        for header, col in required_headers.items():
            if (header is "Completeness" or header is "Contamination") and col is None:
                raise GenomeDatabaseError(
                    "Unable to find %s header in the CheckM file. Check that the CheckM file is correct: %s." % (header, checkm_fh.name))

        return required_headers

    def _processCheckM(self, checkm_file):
        """Parse information from CheckM file.

        Parameters
        ----------
        checkm_file : str
            Name of file containing CheckM results.

        Returns
        -------
        dict
            CheckM statistics for each genome.
        """

        try:
            checkm_fh = open(checkm_file, "rb")
        except:
            raise GenomeDatabaseError(
                "Cannot open CheckM file: " + checkm_file)

        required_headers = self._identifyHeadersCheckM(checkm_fh)

        # populate CheckM results dict
        checkm_results_dict = {}

        for line in checkm_fh:
            line = line.rstrip()
            splitline = line.split("\t")

            bin_id = splitline[required_headers["Bin Id"]]
            completeness = splitline[required_headers["Completeness"]]
            contamination = splitline[required_headers["Contamination"]]
            lineage = splitline[required_headers["Marker lineage"]]
            genome_count = splitline[required_headers["# genomes"]]
            marker_count = splitline[required_headers["# markers"]]
            set_count = splitline[required_headers["# marker sets"]]
            heterogeneity = splitline[required_headers["Strain heterogeneity"]]

            checkm_results_dict[bin_id] = {"completeness": completeness,
                                           "contamination": contamination,
                                           "lineage": lineage,
                                           "genome_count": genome_count,
                                           "marker_count": marker_count,
                                           "set_count": set_count,
                                           "heterogeneity": heterogeneity}

        checkm_fh.close()

        return checkm_results_dict


    # TODO: This should not be here, techincally the backend is agnostic so
    # shouldn't assume command line.
    def _confirm(self, msg):
        raw = raw_input(msg + " (y/N): ")
        if raw.upper() == "Y":
            return True
        return False


    def deleteGenomes(self, batchfile=None, external_ids=None, reason=None):
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

            genome_ids = self._externalGenomeIdsToGenomeIds(external_ids)

            if genome_ids is False:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Unable to retrieve genome ids.")

            has_permission, username, genomes_owners = self._hasPermissionToEditGenomes(
                genome_ids)

            print genome_ids

            if has_permission is None:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Unable to retrieve permissions for genomes.")

            if has_permission is False:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Insufficient permissions.")

            if not self._confirm("Are you sure you want to delete %i genomes (this action cannot be undone)" % len(genome_ids)):
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
    def _externalGenomeIdsToGenomeIds(self, external_ids):
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
    def _hasPermissionToEditGenomes(self, genome_ids):
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
