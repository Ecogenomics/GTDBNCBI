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

from biolib.checksum import sha256

import ConfigMetadata


class TigrfamSearch(object):
    """Runs TIGRfam HMMs over a set of genomes."""

    def __init__(self, cur, currentUser, threads):
        """Initialization."""

        self.cur = cur
        self.currentUser = currentUser

        self.threads = threads
        self.tigrfam_hmms = ConfigMetadata.TIGRFAM_HMMS
        self.protein_file_suffix = ConfigMetadata.PROTEIN_FILE_SUFFIX
        self.tigrfam_suffix = ConfigMetadata.TIGRFAM_SUFFIX
        self.tigrfam_top_hit_suffix = ConfigMetadata.TIGRFAM_TOP_HIT_SUFFIX
        self.checksum_suffix = ConfigMetadata.CHECKSUM_SUFFIX

    def _topHit(self, tigrfam_file):
        """Determine top hits to TIGRFAMs.

        A gene is assigned to a single TIGRFAM
        family. This will be the top hit among
        all TIGRFAM HMMs and pass the threshold
        for the HMM.

        Parameters
        ----------
        tigrfam_file : str
            Name of file containing hits to TIGRFAM HMMs.
        """
        assembly_dir, filename = os.path.split(tigrfam_file)
        output_tophit_file = os.path.join(assembly_dir, filename.replace(self.tigrfam_suffix,
                                                                         self.tigrfam_top_hit_suffix))

        tophits = {}
        for line in open(tigrfam_file):
            if line[0] == '#':
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[3]
            evalue = float(line_split[4])
            bitscore = float(line_split[5])
            if gene_id in tophits:
                if bitscore > tophits[gene_id][2]:
                    tophits[gene_id] = (hmm_id, evalue, bitscore)
            else:
                tophits[gene_id] = (hmm_id, evalue, bitscore)

        fout = open(output_tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, stats in tophits.iteritems():
            hit_str = ','.join(map(str, stats))
            fout.write('%s\t%s\n' % (gene_id, hit_str))
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
            if gene_file is None:
                break

            assembly_dir, filename = os.path.split(gene_file)
            output_hit_file = os.path.join(assembly_dir, filename.replace(self.protein_file_suffix,
                                                                          self.tigrfam_suffix))

            hmmsearch_out = os.path.join(assembly_dir, filename.replace(self.protein_file_suffix, '_tigrfam.out'))
            cmd = 'hmmsearch -o %s --tblout %s --noali --notextw --cut_nc --cpu %d %s %s' % (hmmsearch_out,
                                                                                             output_hit_file,
                                                                                             self.cpus_per_genome,
                                                                                             self.tigrfam_hmms,
                                                                                             gene_file)
            os.system(cmd)

            # calculate checksum
            checksum = sha256(output_hit_file)
            fout = open(output_hit_file + self.checksum_suffix, 'w')
            fout.write(checksum)
            fout.close()

            # identify top hit for each gene
            self._topHit(output_hit_file)

            # allow results to be processed or written to file
            queueOut.put(gene_file)

    def _writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a is None:
                break

            processedItems += 1
            statusStr = '==> Finished processing %d of %d (%.2f%%) genomes.' % (processedItems,
                                                                                numDataItems,
                                                                                float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            #===============================================================================
            # During the TigrFam and Pfam searches no call is made to the database.
            # The transaction is in an open state but idle and after a certain time Postgres or a Firewall is stopping the connection (Timeout error).
            # When the searches are done , GTDB try to submit data to Watson thinking the connection is still open and it crashes.
            # This error does not happen for smaller datasets because the searches are way quicker and pass through the Firewall/PostgresQL timeout threshold.
            # To avoid this crash, we are running a dummy request (SELECT 1) every 'n' searches (here 200) to keep the connection open
            # by resetting the timeout connection countdown to 0.
            #===============================================================================
            if (processedItems % 200) == 0:
                try:
                    self.cur.execute("SELECT 1")
                except:
                    raise

        self.cur.execute("SELECT 1")
        sys.stdout.write('\n')

    def run(self, gene_files):
        """Annotate genes with TIGRFAM HMMs.

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
