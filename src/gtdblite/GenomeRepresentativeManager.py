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

import psycopg2

from biolib.parallel import Parallel

from MarkerSetManager import MarkerSetManager
from AlignedMarkerManager import AlignedMarkerManager


class GenomeRepresentativeManager(object):
    ''''Manage genome representatives.'''

    def __init__(self, cur, currentUser, threads):
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

        # threshold used to assign genome to representative
        self.repThreshold = 0.99

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
            Required AAI to cluster sequences.

        Returns
        -------
        bool
            True if AAI is greater than or equal to threshold, else False.
        """

        assert len(seq1) == len(seq2)

        max_mismatches = (1.0 - threshold) * len(seq1)

        mismatches = 0
        matches = 0
        for c1, c2 in itertools.izip(seq1, seq2):
            if c1 == '-' or c2 == '-':
                continue
            elif c1 != c2:
                mismatches += 1
                if mismatches >= max_mismatches:
                    return False
            else:
                matches += 1

        aai = float(matches) / max(1, (matches + mismatches))

        return aai >= threshold

    def _unprocessedGenomes(self):
        """Identify genomes that have not been compared to representatives.

        Returns
        -------
        list
            List of database identifiers for unprocessed genomes.
        """

        self.cur.execute("SELECT id " +
                         "FROM metadata_taxonomy " +
                         "WHERE gtdb_representative IS NULL")
        unprocessed_genome_ids = [genome_id[0] for genome_id in self.cur.fetchall()]

        return unprocessed_genome_ids

    def representativeGenomes(self):
        """Get list of representative genomes.

        Returns
        -------
        list
            List of database identifiers for representative genomes.
        """

        self.cur.execute("SELECT id " +
                         "FROM metadata_taxonomy " +
                         "WHERE gtdb_representative = 'TRUE'")
        rep_genome_ids = [genome_id[0] for genome_id in self.cur.fetchall()]

        return rep_genome_ids

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

        # identify and align genes from canonical bacterial and archaeal marker sets
        all_markers = set(bac_marker_ids).union(ar_marker_ids)
        aligned_mngr = AlignedMarkerManager(self.cur, self.threads)
        aligned_mngr.calculateAlignedMarkerSets(unprocessed_genome_ids, all_markers)

        # get list of representative genomes
        rep_genome_ids = self.representativeGenomes()
        self.logger.info("Comparing %d unprocessed genomes to %d representatives." %
                                                    (len(unprocessed_genome_ids),
                                                     len(rep_genome_ids)))

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
            rep_bac_aligns[rep_id] = marker_set_mngr.concatenatedAlignedMarkers(rep_id, bac_marker_index)
            rep_ar_aligns[rep_id] = marker_set_mngr.concatenatedAlignedMarkers(rep_id, ar_marker_index)

        # process each genome
        assigned_to_rep_count = 0
        for genome_id in unprocessed_genome_ids:
            # get canonical alignment
            genome_bac_align = marker_set_mngr.concatenatedAlignedMarkers(genome_id, bac_marker_index)
            genome_ar_align = marker_set_mngr.concatenatedAlignedMarkers(genome_id, ar_marker_index)

            for rep_id in rep_genome_ids:
                rep_bac_align = rep_bac_aligns[rep_id]
                rep_ar_align = rep_ar_aligns[rep_id]

                bCluster = (self._aai_test(genome_bac_align, rep_bac_align, self.repThreshold) or
                            self._aai_test(genome_ar_align, rep_ar_align, self.repThreshold))

                if bCluster:
                    # assign genome to current representative
                    assigned_to_rep_count += 1
                    query = ("UPDATE metadata_taxonomy " +
                             "SET gtdb_genome_representative = %s " +
                             "WHERE id = %s")
                    self.cur.execute(query, (rep_id, genome_id))
                    break

        self.logger.info("Assigned %d genomes to a representative." % assigned_to_rep_count)

        # currently, new genomes are never made a representative
        query = "UPDATE metadata_taxonomy SET gtdb_representative = %s WHERE id = %s"
        self.cur.executemany(query, [('False', genome_id) for genome_id in unprocessed_genome_ids])
