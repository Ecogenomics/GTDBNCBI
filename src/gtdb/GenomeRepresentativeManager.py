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
import itertools
import tempfile
import psycopg2
import operator
import os
import shutil
import pickle

from biolib.parallel import Parallel
from biolib.common import remove_extension, make_sure_path_exists

from .Tools import fastaPathGenerator, splitchunks

from .Exceptions import GenomeDatabaseError
from .GenomeManager import GenomeManager
from .MarkerSetManager import MarkerSetManager
from .AlignedMarkerManager import AlignedMarkerManager
from .MetadataManager import MetadataManager

from . import DefaultValues


class GenomeRepresentativeManager(object):
    ''''Manage genome representatives.'''

    def __init__(self, cur, currentUser, threads, db_release):
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
        self.threads = threads

        self.db_release = db_release

        # threshold used to assign genome to representative
        self.aai_threshold = DefaultValues.AAI_CLUSTERING_THRESHOLD
        self.ani_threshold = DefaultValues.FASTANI_CLUSTERING_THRESHOLD

    def _aai_test(self, seq1, seq2, threshold):
        """Test AAI between sequences.

        The identify is only calculate across
        positions where both sequences have
        an amino acid.

        Parameters
        ----------
        seq1 : str
            First sequence.
        seq2 : float
            Second sequence.
        threshold : float
            Minimum AAI required for reporting.

        Returns
        -------
        float
            AAI between sequences if it is greater than the
            specified threshold, else None.
        """

        assert len(seq1) == len(seq2)

        max_mismatches = (1.0 - threshold) * len(seq1)

        mismatches = 0
        matches = 0

        for c1, c2 in zip(seq1, seq2):
            if c1 == '-' or c2 == '-':
                continue
            elif c1 != c2:
                mismatches += 1
                if mismatches >= max_mismatches:
                    return None
            else:
                matches += 1

        aai = float(matches) / max(1, (matches + mismatches))
        if aai < threshold:
            return None

        return aai

    def _aai_mismatches(self, seq1, seq2, max_mismatches):
        """Calculate mismatches between sequences.

        Mismatches are only calculate across
        positions where both sequences have
        an amino acid. The max_mismatches
        threshold is used to quickly stop
        comparisons between divergent
        sequences.

        Parameters
        ----------
        seq1 : str
            First sequence.
        seq2 : str
            Second sequence.
        max_mismatches : int
            Maximum allowed mismatches between sequences.

        Returns
        -------
        int
            Mismatches between sequences if it contains < max_mismatches, else None.
        """

        mismatches = 0
        matches = 0

        for c1, c2 in zip(seq1, seq2):
            if c1 == '-' or c2 == '-':
                continue
            elif c1 != c2:
                mismatches += 1
                if mismatches >= max_mismatches:
                    return None
            else:
                matches += 1

        return mismatches

    def _unprocessedGenomes(self):
        """Identify genomes that have not been compared to representatives.

        Returns
        -------
        list
            List of database identifiers for unprocessed genomes.
        """

        try:
            self.cur.execute("SELECT id " +
                             "FROM metadata_taxonomy " +
                             "WHERE gtdb_representative IS NULL")
            unprocessed_genome_ids = [genome_id[0]
                                      for genome_id in self.cur.fetchall()]

        except GenomeDatabaseError as e:
            raise e

        return unprocessed_genome_ids

    def dereplicatedGenomes(self):
        """Get identifiers from the dereplicated set of genome.

        Identifiers are return for all representative genomes
        and all genomes without a representative. This is the
        typical set of genomes that will be used for tree
        inference and many other downstream applications. No
        filtering for genome quality is performed on this set
        of genomes.

        Returns
        -------
        list
            List of database genome identifiers.
        """

        try:
            self.cur.execute("SELECT id " +
                             "FROM metadata_taxonomy " +
                             "WHERE gtdb_representative = 'TRUE' "
                             "OR gtdb_genome_representative IS NULL")
            derep_genome_ids = [genome_id[0] for genome_id in self.cur]

        except GenomeDatabaseError as e:
            raise e

        return derep_genome_ids

    def ncbiDereplicatedGenomes(self, include_user_reps):
        """Get identifiers from the dereplicated set of NCBI genome.

        Identifiers are return for NCBI representative genomes
        and NCBI genomes without a representative.

        Parameters
        ----------
        include_user_reps : bool
            Flag indicating if NCBI genomes assigned to a
            User representative should be returned.

        Returns
        -------
        list
            List of database genome identifiers.
        """

        try:
            genome_mngr = GenomeManager(self.cur, self.currentUser)
            ncbi_genomes_ids = genome_mngr.ncbiGenomeIds()

            self.cur.execute("SELECT id " +
                             "FROM metadata_taxonomy " +
                             "WHERE (gtdb_representative = 'TRUE' " +
                             "OR gtdb_genome_representative IS NULL) " +
                             "AND id = ANY(%s)", (ncbi_genomes_ids,))
            derep_genome_ids = [genome_id[0] for genome_id in self.cur]

            if include_user_reps:
                self.cur.execute("SELECT id " +
                                 "FROM metadata_taxonomy " +
                                 "WHERE (gtdb_representative = 'FALSE' " +
                                 "AND gtdb_genome_representative LIKE %s " +
                                 "AND id = ANY(%s))", ('U%', ncbi_genomes_ids,))
                derep_genome_ids += [genome_id[0] for genome_id in self.cur]

        except GenomeDatabaseError as e:
            raise e

        return derep_genome_ids

    def sraRepresentatives(self):
        """Get identifiers for representatives SRA genomes.
        Returns
        -------
        list
            List of database genome identifiers.
        """

        try:

            self.cur.execute("SELECT id from sra_dereplicated")
            sra_genome_ids = [genome_id[0] for genome_id in self.cur]

        except GenomeDatabaseError as e:
            raise e

        return sra_genome_ids

    def representativeGenomes(self):
        """Get genome identifiers for all representative genomes.

        Returns
        -------
        list
            List of database identifiers for representative genomes.
        """

        try:
            self.cur.execute("SELECT id " +
                             "FROM metadata_taxonomy " +
                             "WHERE gtdb_representative = 'TRUE'")
            rep_genome_ids = [genome_id[0] for genome_id in self.cur]

        except GenomeDatabaseError as e:
            raise e

        return rep_genome_ids

    def ncbiRepresentativeGenomes(self):
        """Get genome identifiers for all NCBI representative genomes.

        Returns
        -------
        list
            List of database identifiers for NCBI representative genomes.
        """

        try:
            genome_mngr = GenomeManager(self.cur, self.currentUser)
            ncbi_genomes_ids = genome_mngr.ncbiGenomeIds()

            self.cur.execute("SELECT id " +
                             "FROM metadata_taxonomy " +
                             "WHERE gtdb_representative = 'TRUE' " +
                             "AND id = ANY(%s)", (ncbi_genomes_ids,))
            rep_genome_ids = [genome_id[0] for genome_id in self.cur]

        except GenomeDatabaseError as e:
            raise e

        return rep_genome_ids

    def userRepresentativeGenomes(self):
        """Get genome identifiers for all user representative genomes.

        Returns
        -------
        list
            List of database identifiers for user representative genomes.
        """

        try:
            genome_mngr = GenomeManager(self.cur, self.currentUser)
            user_genomes_ids = genome_mngr.userGenomeIds()

            self.cur.execute("SELECT id " +
                             "FROM metadata_taxonomy " +
                             "WHERE gtdb_representative = 'TRUE' " +
                             "AND id = ANY(%s)", (user_genomes_ids,))
            user_rep_genome_ids = [genome_id[0] for genome_id in self.cur]

        except GenomeDatabaseError as e:
            raise e

        return user_rep_genome_ids

    def _getRepresentativeDomain(self):
        """Get genome identifiers and gtdb_domain for all representative genomes.

        Returns
        -------
        dict
            Dictionary of database {identifiers:domain} for representative genomes.
        """

        try:
            self.cur.execute("SELECT id,gtdb_domain " +
                             "FROM metadata_taxonomy " +
                             "WHERE gtdb_representative = 'TRUE'")
            rep_genome_dictionary = {genome_id: gtdb_domain for (
                genome_id, gtdb_domain) in self.cur}

        except GenomeDatabaseError as e:
            raise e

        return rep_genome_dictionary

    def _domainAssignment(self, genome_id, len_arc_marker, len_bac_marker):
        """Assign genome to domain based on present/absence of canonical marker genes."""

        query_al_mark = ("SELECT count(*) " +
                         "FROM aligned_markers am " +
                         "LEFT JOIN marker_set_contents msc ON msc.marker_id = am.marker_id " +
                         "WHERE genome_id = %s and msc.set_id = %s and (evalue <> '') IS TRUE;")

        self.cur.execute(query_al_mark, (genome_id, 1))
        aligned_bac_count = self.cur.fetchone()[0]

        self.cur.execute(query_al_mark, (genome_id, 2))
        aligned_arc_count = self.cur.fetchone()[0]

        arc_aa_per = (aligned_arc_count * 100.0 / len_arc_marker)
        bac_aa_per = (aligned_bac_count * 100.0 / len_bac_marker)

        if arc_aa_per < DefaultValues.DEFAULT_DOMAIN_THRESHOLD and bac_aa_per < DefaultValues.DEFAULT_DOMAIN_THRESHOLD:
            domain = None
        elif bac_aa_per >= arc_aa_per:
            domain = "d__Bacteria"
        else:
            domain = "d__Archaea"

        return domain, arc_aa_per, bac_aa_per

    def domainAssignmentReport(self, outfile):
        """Reports results of automated domain assignment."""

        # identify genomes that have not been compared to representatives
        self.cur.execute("SELECT id FROM metadata_taxonomy")
        genome_ids = [genome_id[0] for genome_id in self.cur.fetchall()]

        # get concatenated alignments for all representatives
        self.cur.execute(
            "SELECT count(*) from marker_set_contents where set_id = 1;")
        len_bac_marker = self.cur.fetchone()[0]

        self.cur.execute(
            "SELECT count(*) from marker_set_contents where set_id = 2;")
        len_arc_marker = self.cur.fetchone()[0]

        # get mapping from internal to external genome IDs
        genome_mngr = GenomeManager(self.cur, self.currentUser)
        external_genome_ids = genome_mngr.genomeIdsToExternalGenomeIds(
            genome_ids)

        # process each genome
        fout = open(outfile, 'w')
        fout.write(
            'Genome Id\tPredicted domain\tArchaeal Marker Percentage\tBacterial Marker Percentage\tNCBI taxonomy\tGTDB taxonomy\n')

        query_taxonomy_req = ("SELECT id, ncbi_taxonomy, gtdb_taxonomy " +
                              "FROM metadata_taxonomy LEFT JOIN gtdb_taxonomy_view USING (id);")
        self.cur.execute(query_taxonomy_req)
        for genome_id, ncbi_taxonomy, gtdb_taxonomy in self.cur.fetchall():
            domain, arc_aa_per, bac_aa_per = self._domainAssignment(
                genome_id, len_arc_marker, len_bac_marker)

            external_genome_id = external_genome_ids[genome_id]
            fout.write('%s\t%s\t%.2f\t%.2f\t%s\t%s\n' % (
                external_genome_id, domain, arc_aa_per, bac_aa_per, ncbi_taxonomy, gtdb_taxonomy))

        fout.close()

    def _calculate_fastani_distance(self, user_genome, genome_reps):
        """ Calculate the FastANI distance between all user genomes and the reference to classify them at the species level

        Parameters
        ----------
        user_leaf : User genome
        genome_reps : list of representatives genomes

        """
        try:
            self.tmp_output_dir = tempfile.mkdtemp()
            make_sure_path_exists(self.tmp_output_dir)

            # we write the two input files for fastani, the query file and
            # reference file
            query_list_file = open(os.path.join(
                self.tmp_output_dir, 'query_list.txt'), 'w')

            # We need to rebuild the path for each unprocessed genomes
            genome_dirs_query = ("SELECT g.id, g.fasta_file_location,gs.external_id_prefix "
                                 "FROM genomes g " +
                                 "LEFT JOIN genome_sources gs ON gs.id = g.genome_source_id " +
                                 "WHERE g.id in %s")
            self.cur.execute(genome_dirs_query,
                             (tuple([user_genome]),))
            raw_results = self.cur.fetchall()
            genome_dir_user = {a: fastaPathGenerator(
                b, c) for a, b, c in raw_results}
            for _k, v in list(genome_dir_user.items()):
                query_list_file.write('{}\n'.format(v))
            query_list_file.close()

            # We need to rebuild the path for each potential reps
            genome_dirs_query = ("SELECT g.id, g.fasta_file_location,gs.external_id_prefix "
                                 "FROM genomes g " +
                                 "LEFT JOIN genome_sources gs ON gs.id = g.genome_source_id " +
                                 "WHERE g.id in %s")
            self.cur.execute(genome_dirs_query,
                             (tuple(list(zip(*genome_reps))[0]),))
            raw_results = self.cur.fetchall()
            genome_dirs = {a: fastaPathGenerator(
                b, c) for a, b, c in raw_results}
            ref_list_file = open(os.path.join(
                self.tmp_output_dir, 'ref_list.txt'), 'w')
            for _k, v in list(genome_dirs.items()):
                ref_list_file.write('{}\n'.format(v))
            ref_list_file.close()

            # run fastANI
            if not os.path.isfile(os.path.join(self.tmp_output_dir, 'query_list.txt')) or not os.path.isfile(os.path.join(self.tmp_output_dir, 'ref_list.txt')):
                raise

            cmd = 'fastANI --ql {0} --rl {1} -o {2} > /dev/null 2>{3}'.format(os.path.join(self.tmp_output_dir, 'query_list.txt'),
                                                                              os.path.join(
                                                                                  self.tmp_output_dir, 'ref_list.txt'),
                                                                              os.path.join(
                                                                                  self.tmp_output_dir, 'results.tab'),
                                                                              os.path.join(self.tmp_output_dir, 'error.log'))
            os.system(cmd)

            if not os.path.isfile(os.path.join(self.tmp_output_dir, 'results.tab')):
                errstr = 'FastANI has stopped:\n'
                if os.path.isfile(os.path.join(self.tmp_output_dir, 'error.log')):
                    with open(os.path.join(self.tmp_output_dir, 'error.log')) as debug:
                        for line in debug:
                            finalline = line
                        errstr += finalline
                raise ValueError(errstr)

            dict_parser_distance = self._parse_fastani_results(
                os.path.join(self.tmp_output_dir, 'results.tab'), genome_dirs, user_genome)
            if len(dict_parser_distance) == 0:
                return None
            sorted_dict = sorted(iter(list(dict_parser_distance.get(
                user_genome).items())), key=lambda _x_y: _x_y[1]['ani'], reverse=True)
            fastani_matching_reference = sorted_dict[0][0]
            shutil.rmtree(self.tmp_output_dir)
            return fastani_matching_reference

        except ValueError as error:
            if os.path.exists(self.tmp_output_dir):
                shutil.rmtree(self.tmp_output_dir)
            raise error
        except Exception as error:
            if os.path.exists(self.tmp_output_dir):
                shutil.rmtree(self.tmp_output_dir)
            raise error

    def _parse_fastani_results(self, fastout_file, genome_dirs, unprocessed_genomes):
        """ Parse the fastani output file


        Parameters
        ----------
        fastout_file : fastani output file.


        Returns
        -------
        dictionary
            dict_results[user_g]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        """
        dict_results = {}
        with open(fastout_file) as fastfile:
            for line in fastfile:
                info = line.strip().split()
                path_genome = info[1]
                ref_genome = None
                for k, v in list(genome_dirs.items()):
                    if v == path_genome:
                        ref_genome = k
                        break
                ani = round(float(info[2]) / 100, 5)
                af = round(float(info[3]) / float(info[4]), 2)
                if unprocessed_genomes in dict_results and self.ani_threshold <= ani:
                    dict_results[unprocessed_genomes][ref_genome] = {
                        "ani": ani, 'af': af}
                elif self.ani_threshold <= ani:
                    dict_results[unprocessed_genomes] = {
                        ref_genome: {"ani": ani, "af": af}}

        return dict_results

    def assignToRepresentative(self):
        """Assign genomes to representatives.

        This method assumes any genomes to process
        have already been committed to the database
        as identification and alignment of canonical
        marker genes is done in independent database
        transactions.
        """

        # identify genomes that have not been compared to representatives
        unprocessed_genome_ids = self._unprocessedGenomes()
        if not unprocessed_genome_ids:
            return

        # get canonical bacterial and archaeal markers
        marker_set_mngr = MarkerSetManager(self.cur, self.currentUser)
        bac_marker_ids = marker_set_mngr.canonicalBacterialMarkers()
        ar_marker_ids = marker_set_mngr.canonicalArchaealMarkers()

        # identify and align genes from canonical bacterial and archaeal marker
        # sets
        all_markers = set(bac_marker_ids).union(ar_marker_ids)
        aligned_mngr = AlignedMarkerManager(
            self.cur, self.threads, self.db_release)
        aligned_mngr.calculateAlignedMarkerSets(
            unprocessed_genome_ids, all_markers)

        # get list of representative genomes
        rep_genome_ids = self.representativeGenomes()
        self.logger.info("Comparing %d unprocessed genomes to %d representatives." %
                         (len(unprocessed_genome_ids),
                          len(rep_genome_ids)))

        # get external genome IDs for representative genomes
        genome_mngr = GenomeManager(self.cur, self.currentUser)
        external_ids = genome_mngr.genomeIdsToExternalGenomeIds(rep_genome_ids)

        # get domains for all representatives
        rep_genome_dictionary = self._getRepresentativeDomain()

        # define desired order of marker genes
        # (order doesn't matter, but must be consistent between genomes)
        bac_marker_index = {}
        for i, marker_id in enumerate(bac_marker_ids):
            bac_marker_index[marker_id] = i

        ar_marker_index = {}
        for i, marker_id in enumerate(ar_marker_ids):
            ar_marker_index[marker_id] = i

        # get concatenated alignments for all representatives
        rep_bac_aligns = {}
        rep_ar_aligns = {}
        for rep_id in rep_genome_ids:
            rep_bac_aligns[rep_id] = marker_set_mngr.concatenatedAlignedMarkers(
                rep_id, bac_marker_index)
            rep_ar_aligns[rep_id] = marker_set_mngr.concatenatedAlignedMarkers(
                rep_id, ar_marker_index)

        self.cur.execute(
            "SELECT count(*) from marker_set_contents where set_id = 1;")
        len_bac_marker = self.cur.fetchone()[0]

        self.cur.execute(
            "SELECT count(*) from marker_set_contents where set_id = 2;")
        len_arc_marker = self.cur.fetchone()[0]
        # process each genome
        assigned_to_rep_count = 0
        for genome_id in unprocessed_genome_ids:

            # get canonical alignment
            genome_bac_align = marker_set_mngr.concatenatedAlignedMarkers(
                genome_id, bac_marker_index)
            genome_ar_align = marker_set_mngr.concatenatedAlignedMarkers(
                genome_id, ar_marker_index)

            domain, _arc_aa_per, _bac_aa_per = self._domainAssignment(
                genome_id, len_arc_marker, len_bac_marker)

            assigned_representative = None
            assigned_representative_dic = {}
            bac_max_mismatches = (1.0 - self.aai_threshold) * \
                (len(genome_bac_align) - genome_bac_align.count('-'))
            ar_max_mismatches = (1.0 - self.aai_threshold) * \
                (len(genome_ar_align) - genome_ar_align.count('-'))
            for rep_id in rep_genome_ids:
                rep_bac_align = rep_bac_aligns[rep_id]
                rep_ar_align = rep_ar_aligns[rep_id]
                if rep_genome_dictionary[rep_id] == domain:
                    if rep_genome_dictionary[rep_id] == 'd__Bacteria':
                        m = self._aai_mismatches(
                            genome_bac_align, rep_bac_align, bac_max_mismatches)
                        if m is not None:  # necessary to distinguish None and 0
                            #assigned_representative = rep_id
                            assigned_representative_dic[rep_id] = m
                            #bac_max_mismatches = m

                    elif rep_genome_dictionary[rep_id] == 'd__Archaea':
                        m = self._aai_mismatches(
                            genome_ar_align, rep_ar_align, ar_max_mismatches)
                        if m is not None:  # necessary to distinguish None and 0
                            assigned_representative_dic[rep_id] = m
                            ar_max_mismatches = m

            # assign genome to current representative
            if assigned_representative_dic:
                sorted_reps = sorted(list(assigned_representative_dic.items(
                )), key=operator.itemgetter(1))[0:10]
                try:
                    assigned_representative = self._calculate_fastani_distance(
                        genome_id, sorted_reps)
                except Exception as error:
                    raise GenomeDatabaseError(error.message)
                if assigned_representative:
                    assigned_to_rep_count += 1
                    query = ("UPDATE metadata_taxonomy " +
                             "SET gtdb_genome_representative = %s " +
                             "WHERE id = %s")
                    self.cur.execute(
                        query, (external_ids[assigned_representative], genome_id))
                    query_taxonomy_req = ("SELECT gtdb_class, gtdb_species," +
                                          "gtdb_phylum, gtdb_family, gtdb_domain, gtdb_order, gtdb_genus " +
                                          "FROM metadata_taxonomy WHERE id = %s;")
                    self.cur.execute(query_taxonomy_req, (genome_id,))
                    if all(v is None or v == '' for v in self.cur.fetchone()):
                        query_taxonomy_update = ("UPDATE metadata_taxonomy as mt_newg SET " +
                                                 "gtdb_class = mt_repr.gtdb_class," +
                                                 "gtdb_species = mt_repr.gtdb_species," +
                                                 "gtdb_phylum = mt_repr.gtdb_phylum," +
                                                 "gtdb_family = mt_repr.gtdb_family," +
                                                 "gtdb_domain =  mt_repr.gtdb_domain," +
                                                 "gtdb_order = mt_repr.gtdb_order," +
                                                 "gtdb_genus =  mt_repr.gtdb_genus " +
                                                 "FROM metadata_taxonomy mt_repr " +
                                                 "WHERE mt_repr.id = %s " +
                                                 "AND mt_newg.id = %s")
                        self.cur.execute(query_taxonomy_update,
                                         (assigned_representative, genome_id))
            else:
                query_taxonomy_req = ("SELECT gtdb_class, gtdb_species," +
                                      "gtdb_phylum, gtdb_family, gtdb_domain, gtdb_order, gtdb_genus " +
                                      "FROM metadata_taxonomy WHERE id = %s;")
                self.cur.execute(query_taxonomy_req, (genome_id,))
                if all(v is None or v == '' for v in self.cur.fetchone()):
                    domain, _arc_aa_per, _bac_aa_per = self._domainAssignment(
                        genome_id, len_arc_marker, len_bac_marker)

                    if domain:
                        self.cur.execute("UPDATE metadata_taxonomy " +
                                         "SET gtdb_domain = %s " +
                                         "WHERE id = %s", (domain, genome_id))

        self.logger.info("Assigned %d genomes to a representative." %
                         assigned_to_rep_count)

        # currently, new genomes are never made a representative
        query = "UPDATE metadata_taxonomy SET gtdb_representative = %s WHERE id = %s"
        self.cur.executemany(query, [('False', genome_id)
                                     for genome_id in unprocessed_genome_ids])
