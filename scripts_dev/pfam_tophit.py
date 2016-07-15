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

__prog_name__ = 'pfam_tophit.py'
__prog_desc__ = 'Determine top Pfam hit for each gene.'

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
import argparse
import multiprocessing as mp
from collections import defaultdict

from biolib.checksum import sha256
from biolib.external.execute import check_dependencies


class PfamTopHit(object):
    """Determine top Pfam hit for each gene

    This script assumes that genomes are stored in individual
    directories in the following format:

    <domain>/<genome_id>/<assembly_id>/<assembly_id>_protein.faa

    where <domain> is either 'archaea' or 'bacteria', and there
      may be multiple assembly_id for a given genome_id. These
      typically represent different strains from a species.

    This is the directory structure which results from extract_ncbi.py.
    """

    def __init__(self):
        self.pfam_ext = '_pfam.tsv'

    def __workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            pfam_file = queueIn.get(block=True, timeout=None)
            if pfam_file is None:
                break

            assembly_dir, filename = os.path.split(pfam_file)
            output_tophit_file = os.path.join(assembly_dir, filename.replace(self.pfam_ext, '_pfam_tophit.tsv'))

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
            print output_tophit_file

            # calculate checksum
            checksum = sha256(output_tophit_file)
            fout = open(output_tophit_file + '.sha256', 'w')
            fout.write(checksum)
            fout.close()

            # allow results to be processed or written to file
            queueOut.put(pfam_file)

    def __writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a is None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processedItems, numDataItems, float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, genome_dir, threads):
        # get path to all unprocessed Pfam HMM result files
        print 'Reading Pfam HMM files.'
        pfam_files = []
        for genome_id in os.listdir(genome_dir):
            cur_genome_dir = os.path.join(genome_dir, genome_id)
            if os.path.isdir(cur_genome_dir):
                for assembly_id in os.listdir(cur_genome_dir):
                    assembly_dir = os.path.join(cur_genome_dir, assembly_id)
                    groups = assembly_id.split('_')
                    processed_assembly_id = '_'.join(groups[:2])
                    pfam_tophit_file = os.path.join(assembly_dir, 'prodigal', processed_assembly_id + '_pfam_tophit.tsv')
                    if os.path.exists(pfam_tophit_file):
                        # verify checksum
                        checksum_file = pfam_tophit_file + '.sha256'
                        if os.path.exists(checksum_file):
                            checksum = sha256(pfam_tophit_file)
                            cur_checksum = open(checksum_file).readline().strip()
                            if checksum == cur_checksum:
                                continue

                    pfam_file = os.path.join(assembly_dir, 'prodigal', processed_assembly_id + self.pfam_ext)
                    if os.path.exists(pfam_file):
                        pfam_files.append(pfam_file)

        print '  Number of unprocessed genomes: %d' % len(pfam_files)

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in pfam_files:
            workerQueue.put(f)

        for _ in range(threads):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self.__workerThread, args=(workerQueue, writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target=self.__writerThread, args=(len(pfam_files), writerQueue))

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

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_dir', help='directory containing genomes in individual directories')
    parser.add_argument('-t', '--threads', type=int, help='number of threads', default=1)

    args = parser.parse_args()

    try:
        p = PfamTopHit()
        p.run(args.genome_dir, args.threads)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
