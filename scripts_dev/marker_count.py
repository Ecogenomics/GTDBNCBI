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

__prog_name__ = 'marker_count.py'
__prog_desc__ = 'Create table showing presence and absence of marker.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2017'
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


class MarkerCount(object):
    """Tabulate presence and absence of genes in marker set."""

    def __init__(self):
        pass

    def __workerThread(self, marker_set, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            genome_id, pfam_top_hit_file, tigrfam_tophit_file = queueIn.get(block=True, timeout=None)
            if tigrfam_tophit_file is None:
                break
                
            marker_count = defaultdict(int)
            for marker_file in  [pfam_top_hit_file, tigrfam_tophit_file]:
                with open(marker_file) as f:
                    f.readline()
                    
                    for line in f:
                        line_split = line.strip().split('\t')
                        for hit_info in line_split[1].split(';'):
                            marker_id, evalue, bitscore = hit_info.split(',')
                            
                            if marker_id in marker_set:
                                marker_count[marker_id] += 1

            # write results to file
            queueOut.put((genome_id, marker_count))

    def __writerThread(self, numDataItems, marker_set, output_file, writerQueue):
        """Store or write results of worker threads in a single thread."""
        
        fout = open(output_file, 'w')
        fout.write('Genome ID\t# present\t# single copy')
        for marker_id in marker_set:
            fout.write('\t%s' % marker_id)
        fout.write('\n')
        
        processedItems = 0
        while True:
            genome_id, marker_count = writerQueue.get(block=True, timeout=None)
            if marker_count is None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processedItems, numDataItems, float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            # write out results
            fout.write(genome_id)
            fout.write('\t%d' % sum([1 for v in marker_count.values() if v >= 1]))
            fout.write('\t%d' % sum([1 for v in marker_count.values() if v == 1]))
            for marker_id in marker_set:
                fout.write('\t%d' % marker_count[marker_id])
            fout.write('\n')

        sys.stdout.write('\n')
        
        fout.close()

    def run(self, genome_path_file, marker_file, output_file, threads):
        # get marker set
        marker_set = set()
        for line in open(marker_file):
            marker_set.add(line.split()[0].strip().replace('PFAM_', '').replace('TIGR_', ''))
            
        # get path to all top hit files
        print 'Getting path to top hit files.'
        top_hit_files = []
        for i, line in enumerate(open(genome_path_file)):
            if i % 1000 == 0:
                print 'Processed %d genomes.' % i
                
            gtdb_id, genome_dir = line.strip().split('\t')
            if not os.path.exists(genome_dir):
                # this can be a moving target since User genomes
                # can be removed at any time and NCBI genomes change
                # between releases
                continue 
  
            genome_id = gtdb_id.replace('RS_', '').replace('GB_', '')
            pfam_tophit_file = os.path.join(genome_dir, 'prodigal', genome_id + '_pfam_tophit.tsv')
            tigrfam_tophit_file = os.path.join(genome_dir, 'prodigal', genome_id + '_tigrfam_tophit.tsv')
            
            if not os.path.exists(pfam_tophit_file):
                print 'Missing Pfam top hit file: %s' % pfam_tophit_file
            else:
                top_hit_files.append([gtdb_id, pfam_tophit_file, tigrfam_tophit_file])

        print '  Number of identified top hit files: %d' % len(top_hit_files)

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in top_hit_files:
            workerQueue.put(f)

        for _ in range(threads):
            workerQueue.put((None, None, None))

        try:
            workerProc = [mp.Process(target=self.__workerThread, args=(marker_set, workerQueue, writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target=self.__writerThread, args=(len(top_hit_files), marker_set, output_file, writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put((None, None))
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()

            writeProc.terminate()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_path_file', help='file specifying path to GTDB genomes')
    parser.add_argument('marker_file', help='list of Pfam and TIGRfam IDs comprising marker set')
    parser.add_argument('output_file', help='output file')
    parser.add_argument('-t', '--threads', type=int, help='number of threads', default=1)

    args = parser.parse_args()

    try:
        p = MarkerCount()
        p.run(args.genome_path_file, args.marker_file, args.output_file, args.threads)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
