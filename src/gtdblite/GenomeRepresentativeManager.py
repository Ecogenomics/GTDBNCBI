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

from MarkerSetManager import MarkerSetManager


class GenomeRepresentativeManager(object):
    ''''Manage genome representatives.'''

    def __init__(self, cur, currentUser):
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

    def assignToRepresentative(self, db_genome_ids):
        """Assign genomes to representatives.

        Parameters
        ----------
        db_genome_ids : list
            Unique database identifier of genomes to process.
        """

        # get canonical bacterial and archaeal markers
        marker_set_mngr = MarkerSetManager(self.cur, self.currentUser)
        bac_marker_ids = marker_set_mngr.canonicalBacterialMarkers()
        ar_marker_ids = marker_set_mngr.canonicalArchaealMarkers()

        # get list of representative genomes
        rep_genome_ids = self.representativeGenomes()
        self.logger.info("Comparing genomes to %d representatives." % len(rep_genome_ids))

        # process each genome
        assigned_to_rep_count = 0
        for genome_id in db_genome_ids:
            # get canonical alignment
            genome_bac_align = marker_set_mngr.concatenatedAlignedMarkers(genome_id, bac_marker_ids)
            genome_ar_align = marker_set_mngr.concatenatedAlignedMarkers(genome_id, ar_marker_ids)

            for rep_id in rep_genome_ids:
                rep_bac_align = marker_set_mngr.concatenatedAlignedMarkers(rep_id, bac_marker_ids)
                rep_ar_align = marker_set_mngr.concatenatedAlignedMarkers(rep_id, ar_marker_ids)

                bCluster = (self._aai_test(genome_bac_align, rep_bac_align, self.repThreshold) or
                            self._aai(genome_ar_align, rep_ar_align, self.repThreshold))

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
        self.cur.executemany(query, [('False', genome_id) for genome_id in db_genome_ids])
