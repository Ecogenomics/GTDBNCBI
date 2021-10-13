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

__prog_name__ = 'ani_pairwise.py'
__prog_desc__ = 'Calculate ANI over a set of genomes using ani_calculator'

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
import ntpath
import itertools
import multiprocessing as mp
from collections import defaultdict


class ANI(object):
    def __init__(self):
        pass

    def __workerThread(self, output_dir, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            gf1, gf2 = queueIn.get(block=True, timeout=None)
            if gf1 == None:
                break
                
            gf1_name = ntpath.basename(gf1).replace('_genes.fna', '')
            gf2_name = ntpath.basename(gf2).replace('_genes.fna', '')

            output_file = os.path.join(output_dir, gf1_name + '~' + gf2_name + '.tsv')
            outdir = os.path.join(output_dir, gf1_name + '~' + gf2_name)
            cmd = 'ani_calculator -genome1fna %s -genome2fna %s -outfile %s -outdir %s' % (gf1, gf2, output_file, outdir)
            print(cmd)
            os.system(cmd)

            # allow results to be processed or written to file
            queueOut.put(gf1)

    def __writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a == None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) genome pairs.' % (processedItems, numDataItems, float(processedItems)*100/numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')
        
    def build_table(self, output_dir):
        """Build table with results."""
        
        output_table = os.path.join(output_dir, 'ani_table.tsv')
        
        ani_table = defaultdict(lambda: defaultdict(float))
        of_table = defaultdict(lambda: defaultdict(float))
        gids = set()
        for result_file in os.listdir(output_dir):
            if not result_file.endswith('.tsv'):
                continue
               
            r = os.path.join(output_dir, result_file)
            if r == output_table:
                continue
                
            with open(r) as f:
                f.readline()
                
                for line in f:
                    line_split = line.strip().split('\t')
                    print(r, line_split)
                    
                    gid1 = line_split[0].replace('_genes.fna', '')
                    gid2 = line_split[1].replace('_genes.fna', '')
                    ani12 = float(line_split[2])
                    ani21 = float(line_split[3])
                    of12 = float(line_split[4])
                    of21 = float(line_split[5])
                    
                    print(r, gid1, gid2, ani12, ani21, of12, of21)

                    gids.add(gid1)
                    gids.add(gid2)
                    ani_table[gid1][gid2] = 0.5*(ani12+ani21)
                    ani_table[gid2][gid1] = 0.5*(ani12+ani21)
                    of_table[gid1][gid2] = 0.5*(of12+of21)
                    of_table[gid2][gid1] = 0.5*(of12+of21)
                 
        # write out table
        fout = open(output_table, 'w')
        for gid in sorted(gids):
            fout.write('\t' + gid)
        fout.write('\n')
        
        for r, gid1 in enumerate(sorted(gids)):
            fout.write(gid1)
            for c, gid2 in enumerate(sorted(gids)):
                if r < c:
                    fout.write('\t%.2f' % (100*of_table[gid1][gid2]))
                elif r == c:
                    fout.write('\t-')
                else:
                    fout.write('\t%.2f' % ani_table[gid1][gid2])
                    
                
            fout.write('\n')
        fout.close()

    def run(self, gene_dir, output_dir, threads):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        # read gene files to process
        gene_files = []
        for gene_file in os.listdir(gene_dir):
            if not gene_file.endswith('.fna'):
                continue
                
            gene_files.append(os.path.join(gene_dir, gene_file))
        
        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        num_pairs = 0
        for gf1, gf2 in itertools.combinations(gene_files, 2):
            workerQueue.put((gf1,gf2))
            num_pairs += 1

        for _ in range(threads):
            workerQueue.put((None, None))

        try:
            workerProc = [mp.Process(target = self.__workerThread, args = (output_dir, workerQueue, writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target = self.__writerThread, args = (num_pairs, writerQueue))

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
            
        self.build_table(output_dir)

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gene_dir', help='directory with called genes in nt space')
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('-t', '--threads', type=int, help='number of threads', default=1)

    args = parser.parse_args()

    try:
        p = ANI()
        p.run(args.gene_dir, args.output_dir, args.threads)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
