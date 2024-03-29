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
import logging
import multiprocessing
import tempfile
import subprocess
import shutil
import time
import gzip

import gtdb.ConfigMetadata as ConfigMetadata
from gtdb.GenomeDatabaseConnection import GenomeDatabaseConnection
from gtdb.Tools import fastaPathGenerator, splitchunks

from biolib.seq_io import read_fasta


class AlignedMarkerManager(object):
    ''''Manage the processing of Aligned Markers and querying marker information.'''

    def __init__(self, cur, threads, release):
        """Initialize.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        threads : int
            Number of threads to use for processing.
        """

        self.logger = logging.getLogger()
        self.threads = threads
        self.release = release

        # self.conn = GenomeDatabaseConnection()
        # self.conn.MakePostgresConnection()
        self.cur = cur

        self.tigrfam_suffix = ConfigMetadata.TIGRFAM_SUFFIX
        self.tigrfam_top_hit_suffix = ConfigMetadata.TIGRFAM_TOP_HIT_SUFFIX
        self.pfam_suffix = ConfigMetadata.PFAM_SUFFIX
        self.pfam_top_hit_suffix = ConfigMetadata.PFAM_TOP_HIT_SUFFIX
        self.genome_file_suffix = ConfigMetadata.GENOME_FILE_SUFFIX
        self.protein_file_suffix = ConfigMetadata.PROTEIN_FILE_SUFFIX

    def calculateAlignedMarkerSets(self, db_genome_ids, marker_ids):
        '''
        Run Hmmalign for PFAM and TIGRFAM missing markers

        :param genome_ids: list of genome ids that are used for the tree step
        :param marker_ids: list of marker ids used for the tree building step
        '''

        self.logger.info('Aligning marker genes not already in the database.')
        # return True

        # We need to rebuild the path to each
        genome_dirs_query = ("SELECT g.id, g.genes_file_location,gs.external_id_prefix "
                             "FROM genomes g " +
                             "LEFT JOIN genome_sources gs ON gs.id = g.genome_source_id " +
                             "WHERE g.id in %s")
        self.cur.execute(genome_dirs_query, (tuple(db_genome_ids),))
        raw_results = self.cur.fetchall()
        genome_dirs = {a: fastaPathGenerator(b, c) for a, b, c in raw_results}

        manager = multiprocessing.Manager()
        out_q = manager.Queue()
        procs = []
        nprocs = self.threads
        for item in splitchunks(genome_dirs, nprocs):
            p = multiprocessing.Process(
                target=self._hmmWorker,
                args=(item, marker_ids, out_q))
            procs.append(p)
            p.start()

        # Collect all results into a single result dict. We know how many dicts
        # with results to expect.
        while out_q.empty():
            time.sleep(1)

        # Wait for all worker processes to finish
        for p in procs:
            p.join()

        return True

    def _hmmWorker(self, subdict_genomes, marker_ids, out_q):
        '''
        The worker function, invoked in a process.

        :param subdict_genomes: sub dictionary of genomes
        :param marker_ids: list of marker ids for one or many marker sets
        :param out_q: manager.Queue()

        '''

        for db_genome_id, path in subdict_genomes.items():
            self._runHmmMultiAlign(db_genome_id, path, marker_ids)
        out_q.put("True")
        return True

    def _runHmmMultiAlign(self, db_genome_id, path, marker_ids):
        '''
        Selects markers that are not aligned for a specific genome.

        :param db_genome_id: Selected genome
        :param path: Path to the genomic fasta file for the genome
        :param marker_ids: list of marker ids for the selected sets
        '''

        temp_con = GenomeDatabaseConnection()
        temp_con.MakePostgresConnection(self.release)
        temp_cur = temp_con.cursor()

        # gather information for all marker genes
        final_genome = []
        final_markerid = []
        final_seq = []
        final_multihits = []
        final_evalue = []
        final_bitscore = []
        final_multihits_number = []
        final_multihits_derep = []

        marker_dbs = {"PFAM": self.pfam_top_hit_suffix,
                      "TIGR": self.tigrfam_top_hit_suffix}
        for marker_db, marker_suffix in marker_dbs.iteritems():
            query = ("SELECT m.id_in_database,m.marker_file_location,m.size,m.id " +
                     "FROM genomes as g, markers as m " +
                     "LEFT JOIN marker_databases as md " +
                     "ON md.id=m.marker_database_id " +
                     "WHERE NOT EXISTS (" +
                     "SELECT * FROM aligned_markers as am " +
                     "WHERE am.genome_id = g.id and am.marker_id = m.id) " +
                     "AND g.id = %s " +
                     "AND m.id in %s " +
                     "AND md.external_id_prefix like %s")

            temp_cur.execute(
                query, (db_genome_id, tuple(marker_ids,), marker_db))
            raw_results = temp_cur.fetchall()
            marker_dict_original = {
                a: {"path": b, "size": c, "db_marker_id": d} for a, b, c, d in raw_results}

            # get all gene sequences
            genome_path = str(path)
            tophit_path = genome_path.replace(
                self.protein_file_suffix, marker_suffix)

            # we load the list of all the genes detected in the genome
            protein_file = tophit_path.replace(
                marker_suffix, self.protein_file_suffix)
                
            # *** DHP: HACK to handle moving of data to specific sub-directories
            if False:
                if marker_db == 'PFAM':
                    tophit_path = tophit_path.replace('prodigal', 'prodigal/pfam_33.1_lite')
                    tophit_path = tophit_path.replace('_pfam_lite_tophit', '_pfam_33.1_lite_tophit')
                else:
                    tophit_path = tophit_path.replace('prodigal', 'prodigal/tigrfam_15.0_lite')
                    tophit_path = tophit_path.replace('_tigrfam_lite_tophit.tsv', '_tigrfam_15.0_lite_tophit.tsv')

            all_genes_dict = read_fasta(protein_file, False)

            # Prodigal adds an asterisks at the end of each called genes,
            # These asterisks sometimes appear in the MSA, which can be an
            # issue for some softwares downstream
            for seq_id, seq in all_genes_dict.iteritems():
                if seq[-1] == '*':
                    all_genes_dict[seq_id] = seq[:-1]

            # we store the tophit file line by line and store the
            # information in a dictionary
            open_file = open
            if tophit_path.endswith('.gz'):
                open_file = gzip.open

            with open_file(tophit_path,'rt') as tp:
                # first line is header line
                tp.readline()
                gene_dict = {}
                for line_tp in tp:
                    linelist = line_tp.split("\t")
                    genename = linelist[0]
                    sublist = linelist[1]
                    if ";" in sublist:
                        diff_markers = sublist.split(";")
                    else:
                        diff_markers = [sublist]

                    for each_gene in diff_markers:
                        sublist = each_gene.split(",")
                        markerid = sublist[0]

                        if markerid not in marker_dict_original:
                            continue

                        evalue = sublist[1]
                        bitscore = float(sublist[2].strip())

                        if markerid in gene_dict:

                            oldbitscore = gene_dict.get(
                                markerid).get("bitscore")
                            if oldbitscore < bitscore:
                                list_multihits = gene_dict.get(
                                    markerid).get('multihits_number')
                                list_multihits.append(
                                    all_genes_dict.get(genename))
                                gene_dict[markerid] = {"marker_path": marker_dict_original.get(markerid).get("path"),
                                                       "gene": genename,
                                                       "gene_seq": all_genes_dict.get(genename),
                                                       "evalue": evalue,
                                                       "bitscore": bitscore,
                                                       "multihit": True,
                                                       "db_marker_id": marker_dict_original.get(markerid).get("db_marker_id"),
                                                       "multihits_number": list_multihits}

                            else:
                                gene_dict.get(markerid)["multihit"] = True
                                gene_dict.get(markerid).get('multihits_number').append(
                                    all_genes_dict.get(genename))
                        else:

                            gene_dict[markerid] = {"marker_path": marker_dict_original.get(markerid).get("path"),
                                                   "gene": genename,
                                                   "gene_seq": all_genes_dict.get(genename),
                                                   "evalue": evalue,
                                                   "bitscore": bitscore,
                                                   "multihit": False,
                                                   "db_marker_id": marker_dict_original.get(markerid).get("db_marker_id"),
                                                   "multihits_number": [all_genes_dict.get(genename)]}
            for mid, info in marker_dict_original.iteritems():
                if mid not in gene_dict:
                    final_genome.append(db_genome_id)
                    final_markerid.append(info.get("db_marker_id"))
                    final_seq.append("-" * info.get("size"))
                    final_multihits.append(False)
                    final_evalue.append(None)
                    final_bitscore.append(None)
                    final_multihits_number.append(None)
                    final_multihits_derep.append(None)

            result_aligns = self._runHmmAlign(gene_dict, db_genome_id)
            for result_align in result_aligns:
                final_genome.append(result_align[0])
                final_markerid.append(result_align[1])
                final_seq.append(result_align[2])
                final_multihits.append(result_align[3])
                final_evalue.append(result_align[4])
                final_bitscore.append(result_align[5])
                final_multihits_number.append(int(result_align[6]))
                final_multihits_derep.append(int(result_align[7]))

        if final_genome:
            query = "SELECT upsert_aligned_markers(%s,%s,%s,%s,%s,%s,%s::int[],%s::int[])"
            temp_cur.execute(query, (final_genome,
                                     final_markerid,
                                     final_seq,
                                     final_multihits,
                                     final_evalue,
                                     final_bitscore,
                                     final_multihits_number,
                                     final_multihits_derep))
        temp_con.commit()
        temp_cur.close()
        temp_con.ClosePostgresConnection()

        return True

    def _runHmmAlign(self, marker_dict, genome):
        '''
        Run hmmalign for a set of genes for a specific genome. This is run in a temp folder.

        :param marker_dict: list of markers that need to be aligned
        :param genome: specific genome id

        Returns
        --------------
        List of tuple to be inserted in aligned_markers table
        '''
        result_genomes_dict = []
        hmmalign_dir = tempfile.mkdtemp()
        input_count = 0
        for _markerid, marker_info in marker_dict.iteritems():

            hmmalign_gene_input = os.path.join(
                hmmalign_dir, "input_gene{0}.fa".format(input_count))
            input_count += 1
            out_fh = open(hmmalign_gene_input, 'wb')
            out_fh.write(">{0}\n".format(marker_info.get("gene")))
            out_fh.write("{0}".format(marker_info.get("gene_seq")))
            out_fh.close()
            proc = subprocess.Popen(["hmmalign", "--outformat", "Pfam", marker_info.get(
                "marker_path"), hmmalign_gene_input], stdout=subprocess.PIPE, stderr=subprocess.PIPE,bufsize=1)
            outs=[]
            with proc.stdout:
                for line in iter(proc.stdout.readline, b''):
                    outs.append(line)
            returncode = proc.wait()



            for line in proc.stderr:
                print(line)

            result = self._getAlignedMarker(
                marker_info.get("gene"), outs)

            if len(result) < 1:
                return "TODO"

            result_genomes_dict.append((genome, marker_info.get("db_marker_id"), result, marker_info.get(
                "multihit"), marker_info.get("evalue"), str(marker_info.get("bitscore")), len(marker_info.get('multihits_number')), len(set(marker_info.get('multihits_number')))))
            input_count += 1

        shutil.rmtree(hmmalign_dir)
        return result_genomes_dict

    def _getAlignedMarker(self, hit_name, result_file):
        '''
        Parse the output of Hmmalign

        :param hit_name: gene name
        :param result_file: output file from Hmmalign
        '''
        hit_seq = None
        mask_seq = None

        for line in result_file:
            splitline = line.split(" ", 1)
            if splitline[0] == hit_name.split(" ", 1)[0]:
                rsplitline = line.rsplit(" ", 1)
                hit_seq = rsplitline[-1]
                continue
            if line[0:len("#=GC RF")] == "#=GC RF":
                rsplitline = line.rsplit(" ", 1)
                mask_seq = rsplitline[-1]

        if mask_seq is None:
            raise Exception("Unable to get mask from hmm align result file")

        if hit_seq is None:
            return None

        aligned_marker = ""
        for pos in xrange(0, len(mask_seq)):
            if mask_seq[pos] != 'x':
                continue
            aligned_marker += hit_seq[pos]
        return aligned_marker
