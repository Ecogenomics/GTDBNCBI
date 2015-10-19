#!/usr/bin/env python

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

__prog_name__ = 'pfam_search.py'
__prog_desc__ = 'Run the pfam_search.pl script over a set of genomes.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import multiprocessing as mp

from biolib.checksum import sha256
from biolib.external.execute import check_dependencies


class PfamSearch(object):
    """Runs pfam_search.pl over a set of genomes.

    This is a modified version of the script developped by Donovan
    Original location:
       /srv/whitlam/bio/db/genome_db/scripts/pfam_search.py

    This script takes a list of proteins paths located in the genome directory for users
    """

    def __init__(self):
        check_dependencies(['hmmsearch'])

        self.pfam_hmm_dir = '/srv/db/pfam/27/'
        self.protein_file_ext = '_proteins.faa'

    def __workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            gene_file = queueIn.get(block=True, timeout=None)
            if gene_file == None:
                break

            genome_dir, filename = os.path.split(gene_file)
            output_hit_file = os.path.join(genome_dir, filename.replace(self.protein_file_ext, '_pfam.tsv'))
            cmd = 'pfam_search.pl -outfile %s -cpu 1 -fasta %s -dir %s' % (output_hit_file, gene_file, self.pfam_hmm_dir)
            os.system(cmd)

            # calculate checksum
            checksum = sha256(output_hit_file)
            fout = open(output_hit_file + '.sha256', 'w')
            fout.write(checksum)
            fout.close()

            queueOut.put(gene_file)

    def __writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a == None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processedItems, numDataItems, float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, genes_path_files, threads):
        # get path to all protein files for a a user submission
        # the sha256 test should not be compulsory because each submission is unique
        print 'Reading genomes.'
        gene_files = []
        for gene_path in genes_path_files:
            genome_dir, filename = os.path.split(gene_path)
            pfam_file = os.path.join(genome_dir, filename.replace(self.protein_file_ext, '_pfam.tsv'))
            if os.path.exists(pfam_file):
                # verify checksum
                checksum_file = pfam_file + '.sha256'
                if os.path.exists(checksum_file):
                    checksum = sha256(pfam_file)
                    cur_checksum = open(checksum_file).readline().strip()
                    if checksum == cur_checksum:
                        continue
                    else:
                        os.remove(checksum_file)
                        os.remove(pfam_file)

            if os.path.exists(gene_path):
                gene_files.append(gene_path)

        print '  Number of unprocessed genomes: %d' % len(gene_files)

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in gene_files:
            workerQueue.put(f)

        for _ in range(threads):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self.__workerThread, args=(workerQueue, writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target=self.__writerThread, args=(len(gene_files), writerQueue))

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
