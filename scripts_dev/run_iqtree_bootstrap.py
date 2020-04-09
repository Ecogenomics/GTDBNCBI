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

__prog_name__ = 'run_iqtree.py'
__prog_desc__ = 'Run IQ-TREE in parallel'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2020'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import time
import argparse
import multiprocessing as mp

class IQTREE(object):
    def __init__(self):
        pass

    def __workerThread(self, queue_in, queue_out):
        """Process each data item in parallel."""
        while True:
            idx = queue_in.get(block=True, timeout=None)
            if idx == None:
                break
                
            cmd = 'iqtree -nt {} -s gtdb_r95_ar_concatenated.faa -m LG+C10+F+G -ft gtdb_r95_ar_phylogeny.wag_gamma.tree -pre gtdb_r95_ar122_iqtree.rep{} -bo 1 -safe'.format(
                    self.cpus_per_tree,
                    idx)

            os.system(cmd)

            # allow results to be processed or written to file
            queue_out.put(idx)

    def __writerThread(self, num_reps, writer_queue):
        """Store or write results of worker threads in a single thread."""
        processed = 0
        while True:
            a = writer_queue.get(block=True, timeout=None)
            if a == None:
                break

            processed += 1
            status = 'Finished processing {} of {} ({:.1f}%) replicates.'.format(
                            processed, 
                            num_reps, 
                            float(processed)*100/num_reps)
            sys.stdout.write('{}\r'.format(status))
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self):
        """Run IQ-TREE in parallel."""
        
        start = time.time()
        
        # R95: estimated to by 60GB per tree
        num_reps = 100
        self.num_trees_in_parallel = 10
        self.cpus_per_tree = 9
    
        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for idx in xrange(0, num_reps):
            worker_queue.put(idx)

        for _ in range(self.num_trees_in_parallel):
            worker_queue.put(None)

        try:
            workerProc = [mp.Process(target = self.__workerThread, args = (worker_queue, writer_queue)) for _ in range(self.num_trees_in_parallel)]
            writeProc = mp.Process(target = self.__writerThread, args = (num_reps, writer_queue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writer_queue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()

            writeProc.terminate()
            
            
        end = time.time()
        
        print('Elapsed time: {:.2f} hours'.format((end - start)/(60*60)))

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    try:
        p = IQTREE()
        p.run()
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
