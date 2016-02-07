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
import multiprocessing as mp
from collections import defaultdict

from biolib.checksum import sha256

import ConfigMetadata


class PfamSearch(object):
    """Runs pfam_search.pl over a set of genomes."""

    def __init__(self, threads):
        """Initialization."""

        self.threads = threads

        self.pfam_hmm_dir = ConfigMetadata.PFAM_HMM_DIR
        self.protein_file_suffix = ConfigMetadata.PROTEIN_FILE_SUFFIX
        self.pfam_suffix = ConfigMetadata.PFAM_SUFFIX
        self.pfam_top_hit_suffix = ConfigMetadata.PFAM_TOP_HIT_SUFFIX
        self.checksum_suffix = ConfigMetadata.CHECKSUM_SUFFIX

    def _topHit(self, pfam_file):
        """Determine top hits to PFAMs.

        A gene may be assigned to multiple
        PFAM families from the same clan. The
        search_pfam.pl script takes care of
        most of these issues and here the results
        are simply parsed.

        Parameters
        ----------
        tigrfam_file : str
            Name of file containing hits to TIGRFAM HMMs.
        """

        assembly_dir, filename = os.path.split(pfam_file)
        output_tophit_file = os.path.join(assembly_dir, filename.replace(self.pfam_suffix,
                                                                         self.pfam_top_hit_suffix))

        tophits = defaultdict(dict)
        for line in open(pfam_file):
            if line[0] == '#' or not line.strip():
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[5]
            evalue = float(line_split[12])
            bitscore = float(line_split[11])
            if gene_id in tophits:
                if hmm_id in tophits[gene_id]:
                    if bitscore > tophits[gene_id][hmm_id][1]:
                        tophits[gene_id][hmm_id] = (evalue, bitscore)
                else:
                    tophits[gene_id][hmm_id] = (evalue, bitscore)
            else:
                tophits[gene_id][hmm_id] = (evalue, bitscore)

        fout = open(output_tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, hits in tophits.iteritems():
            hit_str = []
            for hmm_id, stats in hits.iteritems():
                hit_str.append(hmm_id + ',' + ','.join(map(str, stats)))
            fout.write('%s\t%s\n' % (gene_id, ';'.join(hit_str)))
        fout.close()

        # calculate checksum
        checksum = sha256(output_tophit_file)
        fout = open(output_tophit_file + self.checksum_suffix, 'w')
        fout.write(checksum)
        fout.close()

    def _workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            gene_file = queueIn.get(block=True, timeout=None)
            if gene_file == None:
                break

            genome_dir, filename = os.path.split(gene_file)
            output_hit_file = os.path.join(genome_dir, filename.replace(self.protein_file_suffix,
                                                                        self.pfam_suffix))

            cmd = 'pfam_search.pl -outfile %s -cpu %d -fasta %s -dir %s' % (output_hit_file,
                                                                            self.cpus_per_genome,
                                                                            gene_file,
                                                                            self.pfam_hmm_dir)
            os.system(cmd)

            # calculate checksum
            checksum = sha256(output_hit_file)
            fout = open(output_hit_file + self.checksum_suffix, 'w')
            fout.write(checksum)
            fout.close()

            # identify top hit for each gene
            self._topHit(output_hit_file)

            queueOut.put(gene_file)

    def _writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a == None:
                break

            processedItems += 1
            statusStr = '==> Finished processing %d of %d (%.2f%%) genomes.' % (processedItems,
                                                                                numDataItems,
                                                                                float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, gene_files):
        """Annotate genes with Pfam HMMs.

        Parameters
        ----------
        gene_files : iterable
            Gene files in FASTA format to process.
        """

        self.cpus_per_genome = max(1, self.threads / len(gene_files))

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in gene_files:
            workerQueue.put(f)

        for _ in range(self.threads):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self._workerThread, args=(workerQueue, writerQueue)) for _ in range(self.threads)]
            writeProc = mp.Process(target=self._writerThread, args=(len(gene_files), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()

            writeProc.terminate()
            raise
