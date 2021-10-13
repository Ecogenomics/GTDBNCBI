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

__prog_name__ = 'trnascan.py'
__prog_desc__ = 'Run tRNAscan-SE over a set of genomes.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2017'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import ntpath
import argparse
import subprocess
import multiprocessing as mp

from biolib.common import remove_extension
from biolib.external.execute import check_dependencies
from biolib.checksum import sha256


class tRNAScan(object):
    """Runs Run_tRNAScan-SE over a set of genomes."""

    def __init__(self, gbk_arc_assembly_file, gbk_bac_assembly_file, rfq_arc_assembly_file, rfq_bac_assembly_file):
        check_dependencies(['tRNAscan-SE'])

        self.genome_file_ext = '_genomic.fna'

        self.domain_dict = self.parseAssemblySummary(
            gbk_arc_assembly_file, gbk_bac_assembly_file, rfq_arc_assembly_file, rfq_bac_assembly_file)

    def parseAssemblySummary(self, gbk_arc_assembly_file, gbk_bac_assembly_file, rfq_arc_assembly_file, rfq_bac_assembly_file):
        results = {}
        for arcfile in (gbk_arc_assembly_file, rfq_arc_assembly_file):
            with open(arcfile) as arcf:
                arcf.readline()
                for line in arcf:
                    gid = line.split('\t')[0]
                    results[gid] = 'Archaea'
        for bacfile in (gbk_bac_assembly_file, rfq_bac_assembly_file):
            with open(bacfile) as bacf:
                bacf.readline()
                for line in bacf:
                    gid = line.split('\t')[0]
                    results[gid] = 'Bacteria'
        return results

    def __workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            genome_file = queueIn.get(block=True, timeout=None)
            if genome_file == None:
                break

            assembly_dir, filename = os.path.split(genome_file)
            trna_dir = os.path.join(assembly_dir, 'trna')
            genome_id = filename[0:filename.find('_', 4)]

            if not os.path.exists(trna_dir):
                os.makedirs(trna_dir)

            output_file = os.path.join(trna_dir, genome_id + '_trna.tsv')
            log_file = os.path.join(trna_dir, genome_id + '_trna.log')
            stats_file = os.path.join(trna_dir, genome_id + '_trna_stats.tsv')

            domain_flag = '-B'
            if self.domain_dict.get(genome_id) == 'Archaea':
                domain_flag = '-A'

            #cmd = 'tRNAscan-SE %s -q -Q -o %s -m %s -l %s %s' % (domain_flag, output_file, stats_file, log_file, genome_file)
            # os.system(cmd)

            cmd_to_run = ['tRNAscan-SE', domain_flag, '-q', '-Q', '-o',
                          output_file, '-m', stats_file, '-l', log_file, genome_file]
            proc = subprocess.Popen(
                cmd_to_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            # print proc.returncode
            if proc.returncode != 0:
                raise RuntimeError("%r failed, status code %s stdout %r stderr %r" % (
                    cmd_to_run, proc.returncode, stdout, stderr))
            checksum_file = open(output_file + '.sha256', 'w')
            checksum_file.write('{}\n'.format(sha256(output_file)))
            checksum_file.close()

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
                                                                          float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, input_dir, threads):

        genomes_to_consider = None

        # get path to all genome files
        print('Reading genomes.')
        genome_files = []
        for first_three in os.listdir(input_dir):
            onethird_species_dir = os.path.join(input_dir, first_three)
            print(onethird_species_dir)
            if os.path.isfile(onethird_species_dir):
                continue
            for second_three in os.listdir(onethird_species_dir):
                twothird_species_dir = os.path.join(
                    onethird_species_dir, second_three)
                # print twothird_species_dir
                if os.path.isfile(twothird_species_dir):
                    continue
                for third_three in os.listdir(twothird_species_dir):
                    threethird_species_dir = os.path.join(
                        twothird_species_dir, third_three)
                    # print threethird_species_dir
                    if os.path.isfile(threethird_species_dir):
                        continue
                    for complete_name in os.listdir(threethird_species_dir):
                        full_path = os.path.join(
                            threethird_species_dir, complete_name)
                        if os.path.isfile(full_path):
                            continue
                    trna_dir = os.path.join(full_path, 'trna')
                    genome_id = complete_name[0:complete_name.find('_', 4)]

                    trna_file = os.path.join(trna_dir, genome_id + '_trna.tsv')
                    if os.path.exists(trna_file):
                        # verify checksum
                        checksum_file = trna_file + '.sha256'
                        if os.path.exists(checksum_file):
                            checksum = sha256(trna_file)
                            cur_checksum = open(
                                checksum_file).readline().strip()
                            if checksum == cur_checksum:
                                if genomes_to_consider and genome_id in genomes_to_consider:
                                    print('[WARNING] Genome %s is marked as new or modified, but already has tRNAs called.' % genome_id)
                                    print('[WARNING] Genome is being skipped!')
                                continue

                        print('[WARNING] Genome %s has tRNAs called, but an invalid checksum and was not marked for reannotation.' % genome_id)
                        print('[WARNING] Genome will be reannotated.')

                    elif genomes_to_consider and (genome_id not in genomes_to_consider):
                        print('[WARNING] Genome %s has no Pfam annotations, but is also not marked for processing?' % genome_id)
                        print('[WARNING] Genome will be reannotated!')

                    genome_file = os.path.join(
                        full_path, complete_name + self.genome_file_ext)
                    if os.path.exists(genome_file):
                        if os.stat(genome_file).st_size == 0:
                            print('[Warning] Genome file appears to be empty: %s' % genome_file)
                        else:
                            genome_files.append(genome_file)
        print('  Number of unprocessed genomes: %d' % len(genome_files))

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in genome_files:
            workerQueue.put(f)

        for _ in range(threads):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self.__workerThread,
                                     args=(workerQueue, writerQueue))
                          for _ in range(threads)]
            writeProc = mp.Process(target=self.__writerThread,
                                   args=(len(genome_files), writerQueue))

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
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'genome_dir', help='directory containing genomes in individual directories')
    parser.add_argument('--gbk_bac_assembly_file',
                        help='Assembly summary file downloaded from NCBI')
    parser.add_argument('--gbk_arc_assembly_file',
                        help='Assembly summary file downloaded from NCBI')
    parser.add_argument('--rfq_arc_assembly_file',
                        help='Assembly summary file downloaded from NCBI')
    parser.add_argument('--rfq_bac_assembly_file',
                        help='Assembly summary file downloaded from NCBI')
    parser.add_argument('-t', '--threads',
                        help='number of CPUs to use', type=int, default=32)

    args = parser.parse_args()

    try:
        p = tRNAScan(args.gbk_arc_assembly_file, args.gbk_bac_assembly_file,
                     args.rfq_arc_assembly_file, args.rfq_bac_assembly_file)
        p.run(args.genome_dir, args.threads)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
