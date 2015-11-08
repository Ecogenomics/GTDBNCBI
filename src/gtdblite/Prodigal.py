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
import sys
import logging
import shutil
import multiprocessing as mp

from gtdblite import ConfigMetadata

from biolib.external.prodigal import (Prodigal as BioLibProdigal)
from biolib.checksum import sha256


class Prodigal(object):
    """Perform ab initio gene prediction using Prodigal."""

    def __init__(self, threads):
        """Initialize."""

        self.logger = logging.getLogger()

        self.threads = threads

    def _runProdigal(self, fasta_path):
        """Run Prodigal.

        Parameters
        ----------
        fasta_path : str
            Path to FASTA file to process.
        """

        output_dir, fasta_file = os.path.split(fasta_path)
        genome_id = fasta_file[0:fasta_file.rfind('_')]

        prodigal = BioLibProdigal(1, False)
        summary_stats = prodigal.run([fasta_path], output_dir)
        summary_stats = summary_stats[summary_stats.keys()[0]]

        # rename output files to adhere to GTDB conventions
        aa_gene_file = os.path.join(output_dir, genome_id + ConfigMetadata.PROTEIN_FILE_SUFFIX)
        shutil.move(summary_stats.aa_gene_file, aa_gene_file)

        nt_gene_file = os.path.join(output_dir, genome_id + ConfigMetadata.NT_GENE_FILE_SUFFIX)
        shutil.move(summary_stats.nt_gene_file, nt_gene_file)

        gff_file = os.path.join(output_dir, genome_id + ConfigMetadata.GFF_FILE_SUFFIX)
        shutil.move(summary_stats.gff_file, gff_file)

        # save translation table information
        translation_table_file = os.path.join(output_dir, 'prodigal_translation_table.tsv')
        fout = open(translation_table_file, 'w')
        fout.write('%s\t%d\n' % ('best_translation_table', summary_stats.best_translation_table))
        fout.write('%s\t%.2f\n' % ('coding_density_4', summary_stats.coding_density_4 * 100))
        fout.write('%s\t%.2f\n' % ('coding_density_11', summary_stats.coding_density_11 * 100))
        fout.close()

        checksum = sha256(aa_gene_file)
        fout = open(aa_gene_file + ConfigMetadata.CHECKSUM_SUFFIX, 'w')
        fout.write(checksum)
        fout.close()

        return (aa_gene_file, nt_gene_file, gff_file, translation_table_file)

    def _worker(self, out_dict, worker_queue, writer_queue):
        """This worker function is invoked in a process."""

        while True:
            data = worker_queue.get(block=True, timeout=None)
            if data == None:
                break

            (db_genome_id, file_paths) = data

            file_paths["nt_gene_path"] = None
            file_paths["gff_path"] = None
            file_paths["translation_table_path"] = None

            if file_paths.get("aa_gene_path") is None:
                rtn_files = self._runProdigal(file_paths.get("fasta_path"))
                aa_gene_file, nt_gene_file, gff_file, translation_table_file = rtn_files
                file_paths["aa_gene_path"] = aa_gene_file
                file_paths["nt_gene_path"] = nt_gene_file
                file_paths["gff_path"] = gff_file
                file_paths["translation_table_path"] = translation_table_file

            out_dict[db_genome_id] = file_paths
            writer_queue.put(db_genome_id)

    def _writer(self, num_items, writer_queue):
        """Store or write results of worker threads in a single thread."""
        processed_items = 0
        while processed_items < num_items:
            a = writer_queue.get(block=True, timeout=None)
            if a == None:
                break

            processed_items += 1
            statusStr = '==> Finished processing %d of %d (%.2f%%) genomes.' % (processed_items,
                                                                                num_items,
                                                                                float(processed_items) * 100 / num_items)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, genomic_files):
        """Run Prodigal across a set of genomes.

        Parameters
        ----------
        genomic_files : dict
            Dictionary indicating the genomic and gene file for each genome.
        """

        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for db_genome_id, file_paths in genomic_files.iteritems():
            worker_queue.put([db_genome_id, file_paths])

        for _ in range(self.threads):
            worker_queue.put(None)

        try:
            manager = mp.Manager()
            out_dict = manager.dict()

            worker_proc = [mp.Process(target=self._worker, args=(out_dict, worker_queue, writer_queue)) for _ in range(self.threads)]
            writer_proc = mp.Process(target=self._writer, args=(len(genomic_files), writer_queue))

            writer_proc.start()
            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

            writer_queue.put(None)
            writer_proc.join()
        except:
            for p in worker_proc:
                p.terminate()

            writer_proc.terminate()
            raise

        return out_dict
