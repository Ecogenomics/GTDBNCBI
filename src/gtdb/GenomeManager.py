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
import ntpath
import shutil
import logging
import datetime
import sys
import multiprocessing
import psycopg2
import time

import Config
import ConfigMetadata
import Tools

from Exceptions import GenomeDatabaseError
from MetadataManager import MetadataManager
from Prodigal import Prodigal
from TigrfamSearch import TigrfamSearch
from PfamSearch import PfamSearch

from biolib.checksum import sha256
from biolib.common import make_sure_path_exists
from Tools import splitchunks
from psycopg2.extensions import AsIs


class GenomeManager(object):
    """Manage the processing of new genomes and querying genome information."""

    def __init__(self, cur, currentUser, threads=1):
        """Initialize.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        currentUser : User
            Current user of database.
        threads : int
            Number of threads to use for processing.
        """

        self.logger = logging.getLogger()

        self.cur = cur
        self.currentUser = currentUser
        self.threads = threads

        self.defaultGenomeSourceName = 'user'

        self.genomeFileSuffix = ConfigMetadata.GENOME_FILE_SUFFIX
        self.proteinFileSuffix = ConfigMetadata.PROTEIN_FILE_SUFFIX
        self.ntGeneFileSuffix = ConfigMetadata.NT_GENE_FILE_SUFFIX
        self.gffFileSuffix = ConfigMetadata.GFF_FILE_SUFFIX
        self.tranTableFileSuffix = ConfigMetadata.TRANSLATION_TABLE_SUFFIX
        self.checksumSuffix = ConfigMetadata.CHECKSUM_SUFFIX

        self.genomeCopyDir = Config.GTDB_GENOME_USR_DIR
        self.deprecatedUserDir = Config.GTDB_DPRCTD_USR_DIR
        self.deprecatedGBKDir = Config.GTDB_DPRCTD_GBK_DIR
        self.deprecatedRSQDir = Config.GTDB_DPRCTD_RSQ_DIR

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

    def addGenomes(self, checkm_file, batchfile, study_file):
        """Add new genomes to DB.

        Parameters
        ----------
        checkm_file : str
            Name of file containing CheckM results.
        batchfile : str
            Name of file describing genomes to add.
        study_file : str
            Name of file describing study from which genomes were recovered

        Returns
        -------
        list
            List of database genome identifiers of added genomes.
        """

        try:
            self.tmp_output_dir = tempfile.mkdtemp()

            self.logger.info("Parsing Study file.")
            study_id = self._processStudy(study_file)

            self.logger.info("Reading CheckM file.")
            checkm_results_dict = self._processCheckM(checkm_file)

            genomic_files = self._addGenomeBatch(batchfile, self.tmp_output_dir)

            self.logger.info("Running Prodigal to identify genes.")
            prodigal = Prodigal(self.threads)
            file_paths = prodigal.run(genomic_files)

            self.logger.info("Calculating and storing metadata for each genome.")
            manager = multiprocessing.Manager()
            
            progress_queue = multiprocessing.Queue()
            progress_proc = multiprocessing.Process(target = self._progress, args = (len(genomic_files), progress_queue))
            progress_proc.start()
            
            procs = []
            nprocs = self.threads
            for item in splitchunks(genomic_files, nprocs):
                p = multiprocessing.Process(
                    target=self._addGenomesWorker,
                    args=(item, file_paths, checkm_results_dict, study_id, progress_queue))
                procs.append(p)
                p.start()

            # wait for all worker processes to finish
            for p in procs:
                p.join()
                
            progress_queue.put(None)
            progress_proc.join()

            # annotated genes against TIGRfam and Pfam databases
            self.logger.info("Identifying TIGRfam protein families.")
            gene_files = [file_paths[db_genome_id]['aa_gene_path']
                          for db_genome_id in genomic_files]
            tigr_search = TigrfamSearch(self.threads)
            tigr_search.run(gene_files)

            self.logger.info("Identifying Pfam protein families.")
            pfam_search = PfamSearch(self.threads)
            pfam_search.run(gene_files)
        except:
            if os.path.exists(self.tmp_output_dir):
                shutil.rmtree(self.tmp_output_dir)
            raise

        return genomic_files.keys()
        
    def _progress(self, num_genomes, progress_queue):
        """Track progress of large parallel jobs."""
        
        processed_genomes = 0
        while True:
          a = progress_queue.get(block=True, timeout=None)
          if a == None:
            break

          processed_genomes += 1
          statusStr = '==> Finished processing %d of %d (%.2f%%) genomes.' % (processed_genomes, 
                                                                                num_genomes, 
                                                                                float(processed_genomes)*100/num_genomes)
          sys.stdout.write('%s\r' % statusStr)
          sys.stdout.flush()

        sys.stdout.write('\n')

    def _addGenomesWorker(self, genomic_files, file_paths, checkm_results_dict, study_id, progress_queue):
        '''
        The worker function, invoked in a process.

        :param genomic_files: dictionary {genome_id:{checkm_bin_id:value,aa_gene_path:value,fasta_path:value}}
        :param file_paths : dictionary generated from Prodigal
        :param checkm_results_dict: dictionary of checkm results
        :param study_id = study id
        '''
        metadata_mngr = MetadataManager(self.cur, self.currentUser)
        for db_genome_id, values in genomic_files.iteritems():

            self.cur.execute("UPDATE genomes SET study_id = %s WHERE id = %s",
                             (study_id, db_genome_id))

            genome_file_paths = file_paths[db_genome_id]
            output_dir, _file = os.path.split(
                genome_file_paths["aa_gene_path"])

            bin_id = values['checkm_bin_id']
            if bin_id not in checkm_results_dict:
                raise GenomeDatabaseError(
                    "Couldn't find CheckM result for bin %s." % bin_id)

            metadata_mngr.addMetadata(db_genome_id,
                                      genome_file_paths["fasta_path"],
                                      genome_file_paths["gff_path"],
                                      checkm_results_dict[bin_id],
                                      output_dir)
                                      
            progress_queue.put(bin_id)
            
        return True

    def allGenomeIds(self):
        """Get genome identifiers for all genomes.

        Returns
        -------
        list
            Database identifiers for all genomes.
        """

        try:
            query = "SELECT id FROM genomes"
            self.cur.execute(query)

            result_ids = [genome_id[0] for genome_id in self.cur]

        except GenomeDatabaseError as e:
            raise e

        return result_ids

    def userGenomeIds(self):
        """Get genome identifiers for all user genomes.

        Returns
        -------
        list
            Database identifiers for all user genomes.
        """

        try:
            self.cur.execute("SELECT id " +
                             "FROM genome_sources " +
                             "WHERE name = %s", ('user',))

            user_source_id = self.cur.fetchone()[0]

            self.cur.execute("SELECT id " +
                             "FROM genomes " +
                             "WHERE genome_source_id = %s", (user_source_id,))

            result_ids = [genome_id[0] for genome_id in self.cur]

        except GenomeDatabaseError as e:
            raise e

        return result_ids

    def ncbiGenomeIds(self):
        """Get genome identifiers for all NCBI genomes.

        Returns
        -------
        list
            Database identifiers for all NCBI genomes.
        """

        try:
            self.cur.execute("SELECT id " +
                             "FROM genome_sources " +
                             "WHERE name IN %s", (('RefSeq', 'GenBank'),))

            source_ids = [r[0] for r in self.cur.fetchall()]

            self.cur.execute("SELECT id " +
                             "FROM genomes " +
                             "WHERE genome_source_id = ANY(%s)", (source_ids,))

            result_ids = [genome_id[0] for genome_id in self.cur]

        except GenomeDatabaseError as e:
            raise e

        return result_ids

    def moveGenomes(self, db_genome_ids):
        """Move genome files into database directory structure.

        This function assumes addGenomes() has been called. It is
        not directly called by addGenomes() as all database
        queries are performed before moving genomes.

        Parameters
        ----------
        db_genome_ids : list
            Unique database identifiers for genomes.
        """

        assert(self.tmp_output_dir)

        # get database genome identifiers
        self.cur.execute("SELECT genomes.id,user_editable, external_id_prefix || '_' || id_at_source as external_id " +
                         "FROM genomes, genome_sources " +
                         "WHERE genome_source_id = genome_sources.id " +
                         "AND genomes.id in %s", (tuple(db_genome_ids),))

        external_id_dict = {}
        for (genome_id, user_editable, external_id) in self.cur:
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

        gtdb_target_dir = os.path.join(self.genomeCopyDir, username)
        for db_genome_id, external_id in external_id_dict.items():
            tmp_genome_dir = os.path.join(self.tmp_output_dir, external_id)

            genome_target_dir = os.path.join(gtdb_target_dir, external_id)
            if os.path.exists(genome_target_dir):
                raise GenomeDatabaseError(
                    "Genome directory already exists: %s" % genome_target_dir)

            shutil.move(tmp_genome_dir, genome_target_dir)

            self.cur.execute("UPDATE genomes SET fasta_file_location = %s , genes_file_location = %s WHERE id = %s", (
                os.path.join(
                    username, external_id, external_id + self.genomeFileSuffix),
                os.path.join(
                    username, external_id, external_id + self.proteinFileSuffix),
                db_genome_id))

        shutil.rmtree(self.tmp_output_dir)

    def copyGenomes(self, db_genome_ids, genomic, gene, out_dir, gtdb_header):
        """Copy genome data files to specified directory.

        Parameters
        ----------
        db_genome_ids : list
            Unique database identifiers for genomes.
        genomic : boolean
            Flag indicating if genomic data should be copied.
        gene : boolean
            Flag indicating if gene data should be copied.
        out_dir : str
            Output directory for data.
        """

        # get database genome identifiers
        try:
            self.cur.execute("SELECT external_id_prefix, fasta_file_location, genes_file_location " +
                             "FROM genomes, genome_sources " +
                             "WHERE genome_source_id = genome_sources.id " +
                             "AND genomes.id in %s", (tuple(db_genome_ids),))

            for (external_id_prefix, fasta_file_location, genes_file_location) in self.cur:
                dir_prefix = None
                if external_id_prefix == 'U':
                    dir_prefix = Config.GTDB_GENOME_USR_DIR
                elif external_id_prefix == 'RS':
                    dir_prefix = Config.GTDB_GENOME_RSQ_DIR
                elif external_id_prefix == 'GB':
                    dir_prefix = Config.GTDB_GENOME_GBK_DIR
                else:
                    raise GenomeDatabaseError(
                        "Unrecognized database prefix: %s" % external_id_prefix)

                if genomic:
                    genomic_file = os.path.join(dir_prefix, fasta_file_location)
                    if gtdb_header and external_id_prefix != 'U':
                        gtdb_filename = external_id_prefix + "_" + os.path.basename(genomic_file)
                        out_dir_with_header = os.path.join(out_dir, gtdb_filename)
                        shutil.copy(genomic_file, out_dir_with_header)
                    else:
                        shutil.copy(genomic_file, out_dir)

                elif gene:
                    gene_file = os.path.join(dir_prefix, genes_file_location)
                    if gtdb_header and external_id_prefix != 'U':
                        gtdb_filename = external_id_prefix + "_" + os.path.basename(gene_file)
                        out_dir_with_header = os.path.join(out_dir, gtdb_filename)
                        shutil.copy(gene_file, out_dir_with_header)
                    else:
                        shutil.copy(gene_file, out_dir)

        except GenomeDatabaseError as e:
            raise e

    def _addGenomeBatch(self, batchfile, output_dir):
        """Add genomes specific in batch file to DB.

        Parameters
        ----------
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
                self.logger.warning(
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

            abs_fasta_path = os.path.abspath(fasta_path)
            real_fasta_path = os.path.realpath(fasta_path)

            abs_gene_path = None
            if gene_path is not None and gene_path != '':
                abs_gene_path = os.path.realpath(gene_path)

            genome_id = self._addGenomeToDB(
                real_fasta_path, name, desc, source_name, id_at_source, abs_gene_path)
            if not (genome_id):
                raise GenomeDatabaseError(
                    "Failed to add genome: %s" % real_fasta_path)

            self.cur.execute("SELECT external_id_prefix || '_' || id_at_source as external_id " +
                             "FROM genomes, genome_sources " +
                             "WHERE genome_source_id = genome_sources.id " +
                             "AND genomes.id = %s", (genome_id,))

            external_genome_id = self.cur.fetchone()[0]

            genome_output_dir = os.path.join(output_dir, external_genome_id)
            if not os.path.exists(genome_output_dir):
                os.makedirs(genome_output_dir)

            fasta_target_file = os.path.join(
                genome_output_dir, external_genome_id + self.genomeFileSuffix)
            shutil.copy(real_fasta_path, fasta_target_file)

            genes_target_file = None
            if abs_gene_path:
                genes_target_file = os.path.join(
                    output_dir, external_genome_id, external_genome_id + self.proteinFileSuffix)
                shutil.copy(abs_gene_path, genes_target_file)

            genomic_files[genome_id] = {"checkm_bin_id": os.path.splitext(os.path.basename(abs_fasta_path))[0],
                                        "aa_gene_path": genes_target_file,
                                        "fasta_path": fasta_target_file}
        return genomic_files

    def _addGenomeToDB(self, fasta_file_path, name, desc,
                       source, id_at_source, gene_path):
        """Add genome to database.

        Parameters
        ----------
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

            self.cur.execute(
                "SELECT id, external_id_prefix, user_editable FROM genome_sources WHERE name = %s", (source,))
            source_id = None

            for (db_id, _external_id_prefix, user_editable) in self.cur:
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
                # We use update to return a value. This update should fix the concurreny of multit thread using the same value. Update locks the cell during the transaction.
                self.cur.execute("SELECT update_last_auto(%s);", (source_id,))
                id_at_source = str(self.cur.fetchone()[0])

            added = datetime.datetime.now()

            owner_id = None
            if not self.currentUser.isRootUser():
                owner_id = self.currentUser.getUserId()

            self.cur.execute(
                "SELECT id FROM genomes WHERE genome_source_id = %s AND id_at_source = %s", (source_id, id_at_source))

            result = self.cur.fetchall()

            columns = "(name, description, owned_by_root, owner_id, fasta_file_location, " + \
                "fasta_file_sha256, genes_file_location, genes_file_sha256,genome_source_id, id_at_source, date_added)"

            if len(result):
                raise GenomeDatabaseError(
                    "Genome source '%s' already contains id '%s'. Use -f to force an overwrite." % (source, id_at_source))

            self.cur.execute("INSERT INTO genomes " + columns + " "
                             "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) " +
                             "RETURNING id",
                             (name, desc, self.currentUser.isRootUser(), owner_id, fasta_file_path, fasta_sha256_checksum, gene_path, gene_sha256_checksum, source_id, id_at_source, added))
            (db_genome_id,) = self.cur.fetchone()

            return db_genome_id

        except GenomeDatabaseError as e:
            raise e

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

    def _processStudy(self, study_file):
        """Parse information from study file.

        Parameters
        ----------
        study_file : str
            Name of file describing study from which genomes were recovered

        Returns
        -------
        int
            Database identifier of study.
        """

        try:
            study_fh = open(study_file, "rb")
        except:
            raise GenomeDatabaseError(
                "Cannot open study file: " + study_file)

        # required information in study file
        study_info = {
            "study_description": None,
            "sequencing_platform": None,
            "read_files": None,
            "qc_program": None,
            "assembly_program": None,
            "gap_filling_program": None,
            "mapping_program": None,
            "binning_program": None,
            "scaffolding_program": None,
            "genome_assessment_program": None,
            "refinement_description": None,
        }

        for line in study_fh:
            line = line.rstrip()
            if line:
                line_split = line.split("\t")

                if len(line_split) == 1:
                    field = line_split[0]
                    value = ''
                else:
                    field = line_split[0]
                    value = line_split[1]

                if field not in study_info:
                    raise GenomeDatabaseError(
                        "Study file %s contains an unknown field: %s" % (study_file, field))
                else:
                    study_info[field] = value

        study_fh.close()

        # check that all fields were populated
        for field, value in study_info.iteritems():
            if value is None:
                raise GenomeDatabaseError(
                    "Study file %s is missing the field: %s" % (study_file, field))

        # add information to study table
        query = "INSERT INTO study (%s) VALUES %s RETURNING study_id"
        self.cur.execute(
            query, (AsIs(','.join(study_info.keys())), tuple(study_info.values())))
        study_id = self.cur.fetchone()[0]

        return study_id

    # TODO: This should not be here, techincally the backend is agnostic so
    # shouldn't assume command line.
    def _confirm(self, msg):
        raw = raw_input(msg + " (y/N): ")
        if raw.upper() == "Y":
            return True
        return False

    def deleteGenomes(self, batchfile=None, db_genome_ids=None, reason=None):
        '''
        Delete Genomes
        Returns True for success or False for fail

        Parameters:
        :param batchfile: text file listing a range of ids to delete
        :param db_genome_ids: a list of ids can be written directly in the command line
        '''

        self._loggerSetup()

        try:
            if db_genome_ids is False:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Unable to retrieve genome ids.")

            # restrict deletion to genomes owned by user
            has_permission, username, genomes_owners = self._hasPermissionToEditGenomes(
                db_genome_ids)

            if has_permission is None:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Unable to retrieve permissions for genomes.")

            if has_permission is False:
                raise GenomeDatabaseError(
                    "Unable to delete genomes. Insufficient permissions.")

            if db_genome_ids:
                if not self._confirm("Are you sure you want to delete %i genomes (this action cannot be undone)" % len(db_genome_ids)):
                    raise GenomeDatabaseError("User aborted database action.")

                self.cur.execute("DELETE FROM aligned_markers " +
                                 "WHERE genome_id IN %s ", (tuple(db_genome_ids),))

                self.cur.execute("DELETE FROM genome_list_contents " +
                                 "WHERE genome_id IN %s", (tuple(db_genome_ids),))

                # Deletion of metadata

                self.cur.execute("DELETE FROM metadata_genes " +
                                 "WHERE id IN %s", (tuple(db_genome_ids),))
                self.cur.execute("DELETE FROM metadata_ncbi " +
                                 "WHERE id IN %s", (tuple(db_genome_ids),))
                self.cur.execute("DELETE FROM metadata_nucleotide " +
                                 "WHERE id IN %s", (tuple(db_genome_ids),))
                self.cur.execute("DELETE FROM metadata_taxonomy " +
                                 "WHERE id IN %s", (tuple(db_genome_ids),))

                self.cur.execute("DELETE FROM genomes " +
                                 "WHERE id IN %s", (tuple(db_genome_ids),))

                for genome, info in genomes_owners.iteritems():
                    if str(username) != str(info.get("owner")):
                        logging.info('''Genome {0} has been deleted by {1} for the following reason '{2}'
                                          WARNING: {1} is not the owner of this {0} (real owner {3} )
                                          {0} needs to be moved manually to the deprecated folder'''.format(genome, username, reason, info.get("owner")))
                    else:
                        if info.get("prefix") is "U":
                            target = os.path.dirname(
                                os.path.join(self.deprecatedUserDir, info.get("relative_path")))
                        elif info.get("prefix") is "GB":
                            target = os.path.join(
                                self.deprecatedGBKDir, info.get("relative_path"))
                        elif info.get("prefix") is "RS":
                            target = os.path.join(
                                self.deprecatedRSQDir, info.get("relative_path"))
                        make_sure_path_exists(target)
                        os.rename(
                            os.path.dirname(Tools.fastaPathGenerator(info.get("relative_path"), info.get("prefix"))), target)
                        logging.info("Genome {0} has been deleted by {1} for the following reason '{2}'".format(
                            genome, username, reason))
        except GenomeDatabaseError as e:
            raise e

        return True

    def _hasPermissionToEditGenomes(self, db_genome_ids):
        '''
        Function _hasPermissionToEditGenomes
        Check if a user is entitled to delete genomes.
        Users can delete their own genomes, Admin can delete any genomes

        :param db_genome_ids:list of genomes is to delete

        Return a tuple containing:
        - Boolean returning the state of the function
        - The username currently running the delete function
        -a Dictionary listing the list of genomes and where each genome has saved
            the owner,the Prefix (U,GB or RS) , the relative path

        '''
        try:
            if not db_genome_ids:
                raise GenomeDatabaseError(
                    "Unable to retrieve genome permissions, no genomes given: %s" % str(db_genome_ids))

            self.cur.execute("SELECT gs.external_id_prefix,gs.external_id_prefix || '_'|| genomes.id_at_source, owner_id, username, owned_by_root,fasta_file_location "
                             "FROM genomes " +
                             "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
                             "LEFT JOIN genome_sources gs ON gs.id = genomes.genome_source_id " +
                             "WHERE genomes.id in %s", (tuple(db_genome_ids),))

            dict_genomes_user = {}
            for (prefix, public_id, owner_id, username, owned_by_root, fasta_path) in self.cur:

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
        except GenomeDatabaseError as e:
            raise e

        return (True, current_username, dict_genomes_user)

    def genomeIdsToExternalGenomeIds(self, db_genome_ids):
        """Get external genome identifiers from database identifiers.

        Parameters
        ----------
        db_genome_ids : list
            List of database genome ids.

        Returns
        -------
        dict : d[db_genome_id] -> external genome id
            Map from database genome ids to external genome ids.
        """

        external_genome_id = {}
        if db_genome_ids:
            self.cur.execute("SELECT genomes.id, external_id_prefix || '_' || id_at_source as external_id " +
                             "FROM genomes, genome_sources " +
                             "WHERE genome_source_id = genome_sources.id " +
                             "AND genomes.id IN %s", (tuple(db_genome_ids),))

            for internal_id, external_id in self.cur:
                external_genome_id[internal_id] = external_id

        return external_genome_id

    def externalGenomeIdsToGenomeIds(self, external_ids):
        """Get database genome identifiers from external identifiers.

        Parameters
        ----------
        external_ids : list
            List of external genome ids.

        Returns
        -------
        list
            List of database genome ids.
        """

        try:

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
                self.cur.execute("CREATE TEMP TABLE %s (prefix text)" %
                                 (temp_table_name,))
                query = "INSERT INTO {0} (prefix) VALUES (%s)".format(
                    temp_table_name)
                self.cur.executemany(query, [(x,)
                                             for x in map_sources_to_ids.keys()])
            else:
                raise GenomeDatabaseError(
                    "No genome sources found for these ids. %s" % str(external_ids))

            # Find any given tree prefixes that arent in the genome sources
            query = ("SELECT prefix FROM {0} " +
                     "WHERE prefix NOT IN ( " +
                     "SELECT external_id_prefix " +
                     "FROM genome_sources)").format(temp_table_name)

            self.cur.execute(query)

            missing_genome_sources = {}
            for (query_prefix,) in self.cur:
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
                self.cur.execute(
                    "CREATE TEMP TABLE %s (id_at_source text)" % (temp_table_name,))
                query = "INSERT INTO {0} (id_at_source) VALUES (%s)".format(
                    temp_table_name)
                self.cur.executemany(
                    query, [(x,) for x in map_sources_to_ids[source_prefix].keys()])

                # Check to see if there are any that don't exist
                query = ("SELECT id_at_source FROM {0} " +
                         "WHERE id_at_source NOT IN ( " +
                         "SELECT id_at_source " +
                         "FROM genomes, genome_sources " +
                         "WHERE genome_source_id = genome_sources.id " +
                         "AND external_id_prefix = %s)").format(temp_table_name)

                self.cur.execute(query, (source_prefix,))

                missing_ids = []
                for (id_at_source,) in self.cur:
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

                self.cur.execute(query, (source_prefix,))

                for (genome_id,) in self.cur:
                    result_ids.append(genome_id)

        except GenomeDatabaseError as e:
            raise e

        return result_ids

    def printGenomeDetails(self, genome_id_list):
        """Print genome details.

        Parameters
        ----------
        genome_id_list : iterable
            Unique identifier of genomes in database.

        Returns
        -------
        list
            Column headers.
        list
            Content for each row.
        """

        try:
            if not genome_id_list:
                raise GenomeDatabaseError(
                    "Unable to print genomes. No genomes found.")

            columns = "genomes.id, genomes.name, description, owned_by_root, username, " + \
                "external_id_prefix || '_' || id_at_source as external_id, date_added"

            self.cur.execute("SELECT " + columns + " FROM genomes " +
                             "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
                             "JOIN genome_sources AS sources ON genome_source_id = sources.id " +
                             "AND genomes.id in %s " +
                             "ORDER BY genomes.id ASC", (tuple(genome_id_list),))

            header = (
                "genome_id", "name", "description", "owner", "data_added")

            rows = []
            for (_genome_id, name, description, owned_by_root, username, external_id, date_added) in self.cur:
                rows.append((external_id, name, description,
                             ("root" if owned_by_root else username),
                             date_added.date()))

        except GenomeDatabaseError as e:
            raise e

        return header, rows
