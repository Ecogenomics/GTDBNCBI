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

__prog_name__ = 'trnascan_user_genomes.py'
__prog_desc__ = 'Run tRNAscan-SE over User genomes.'

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
import ntpath
import argparse
import multiprocessing as mp

from biolib.common import remove_extension
from biolib.external.execute import check_dependencies
from biolib.external.prodigal import Prodigal
from biolib.checksum import sha256


class tRNAScan(object):
    """Runs Run_tRNAScan-SE over a set of genomes."""

    def __init__(self):
        check_dependencies(['tRNAscan-SE'])
        
        self.genome_file_ext = '_genomic.fna'
        
    def __workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            genome_file, domain = queueIn.get(block=True, timeout=None)
            if genome_file == None:
                break
                
            genome_dir, filename = os.path.split(genome_file)
            trna_dir = os.path.join(genome_dir, 'trna')
            genome_id = filename[0:filename.find('_', 2)]
                
            if os.path.exists(trna_dir):
                queueOut.put(genome_file) # *** Do not preprocess
            else:
                os.makedirs(trna_dir)
                
                output_file = os.path.join(trna_dir, genome_id + '_trna.tsv')
                log_file = os.path.join(trna_dir, genome_id + '_trna.log')
                stats_file = os.path.join(trna_dir, genome_id + '_trna_stats.tsv')
                        
                domain_flag = '-B'
                if domain == 'd__Archaea':
                    domain_flag = '-A'
                
                cmd = 'tRNAscan-SE %s -q -Q -o %s -m %s -l %s %s' % (domain_flag, 
                                                                        output_file, 
                                                                        stats_file, 
                                                                        log_file, 
                                                                        genome_file)
                os.system(cmd)
                
                queueOut.put(genome_file)
        
    def __writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a == None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processedItems, 
                                                                            numDataItems, 
                                                                            float(processedItems)*100/numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')
        
    def _read_domain(self, gtdb_metadata_file):
        """Read domain information from metadata table."""
        
        domain = {}
        with open(gtdb_metadata_file) as f:
            headers = f.readline().strip().split('\t')
            
            domain_index = headers.index('gtdb_domain')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[0]
                domain[gid] = line_split[domain_index]
                
        return domain
            
    def run(self, user_genome_dir, gtdb_metadata_file, threads):
    
        domain = self._read_domain(gtdb_metadata_file)
    
        genomes_to_consider = None

        # get path to all genome files
        print 'Reading genomes.'
        genome_files = []
        for user_id in os.listdir(user_genome_dir):
            cur_user_dir = os.path.join(user_genome_dir, user_id)
            if os.path.isdir(cur_user_dir):
                for genome_id in os.listdir(cur_user_dir):
                    genome_dir = os.path.join(cur_user_dir, genome_id)
                    genome_file = os.path.join(genome_dir, genome_id + self.genome_file_ext)
                    
                    trna_dir = os.path.join(genome_dir, 'trna')
                    trna_file = os.path.join(trna_dir, genome_id + '_trna_stats.tsv')
                    if os.path.exists(trna_file):
                        continue
                    
                    if os.path.exists(genome_file):
                        if os.stat(genome_file).st_size == 0:
                            print '[Warning] Genome file appears to be empty: %s' % genome_file
                        else:
                            if genome_id not in domain:
                                print 'Genome %s is missing domain information: %s' % (genome_id, user_id)
                            else:
                                genome_files.append((genome_file, domain[genome_id]))

        print '  Number of unprocessed genomes: %d' % len(genome_files)

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f, d in genome_files:
            workerQueue.put((f,d))

        for _ in range(threads):
            workerQueue.put((None,None))

        try:
            workerProc = [mp.Process(target = self.__workerThread, 
                            args = (workerQueue, writerQueue)) 
                            for _ in range(threads)]
            writeProc = mp.Process(target = self.__writerThread, 
                                    args = (len(genome_files), writerQueue))

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
  parser.add_argument('user_genome_dir', help='directory containing User genomes')
  parser.add_argument('gtdb_metadata_file', help='GTDB metadata file in TSV format')
  parser.add_argument('-t', '--threads', help='number of CPUs to use', type=int, default=32)

  args = parser.parse_args()

  try:
    p = tRNAScan()
    p.run(args.user_genome_dir, args.gtdb_metadata_file, args.threads)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
