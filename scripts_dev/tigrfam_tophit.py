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

__prog_name__ = 'tigrfam_tophit.py'
__prog_desc__ = 'Determine top TIGRfam hit for each gene.'

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

from biolib.checksum import sha256
from biolib.external.execute import check_dependencies


class TigrfamTopHit(object):
  """Determine top TIGRfam hit for each gene

  This script assumes that genomes are stored in individual
  directories in the following format:

  <domain>/<genome_id>/<assembly_id>/<assembly_id>_protein.faa

  where <domain> is either 'archaea' or 'bacteria', and there
    may be multiple assembly_id for a given genome_id. These
    typically represent different strains from a species.

  This is the directory structure which results from extract_ncbi.py.
  """

  def __init__(self):
    self.tigrfam_ext = '_tigrfam.tsv'

  def __workerThread(self, queueIn, queueOut):
    """Process each data item in parallel."""
    while True:
      tigrfam_file = queueIn.get(block=True, timeout=None)
      if tigrfam_file == None:
        break

      assembly_dir, filename = os.path.split(tigrfam_file)
      output_tophit_file = os.path.join(assembly_dir, filename.replace(self.tigrfam_ext, '_tigrfam_tophit.tsv'))

      tophits = {}
      for line in open(tigrfam_file):
        if line[0] == '#' or line[0] == '[':
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
      fout = open(output_tophit_file + '.sha256', 'w')
      fout.write(checksum)
      fout.close()

      # allow results to be processed or written to file
      queueOut.put(tigrfam_file)

  def __writerThread(self, numDataItems, writerQueue):
    """Store or write results of worker threads in a single thread."""
    processedItems = 0
    while True:
      a = writerQueue.get(block=True, timeout=None)
      if a == None:
        break

      processedItems += 1
      statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processedItems, numDataItems, float(processedItems)*100/numDataItems)
      sys.stdout.write('%s\r' % statusStr)
      sys.stdout.flush()

    sys.stdout.write('\n')

  def run(self, genome_dir, threads):
    # get path to all unprocessed TIGRfam HMM result files
    print 'Reading TIGRfam HMM files.'
    tigrfam_files = []
    for genome_id in os.listdir(genome_dir):
      cur_genome_dir = os.path.join(genome_dir, genome_id)
      if os.path.isdir(cur_genome_dir):
        for assembly_id in os.listdir(cur_genome_dir):
          assembly_dir = os.path.join(cur_genome_dir, assembly_id)

          tigrfam_tophit_file = os.path.join(assembly_dir, assembly_id + '_tigrfam_tophit.tsv')
          if os.path.exists(tigrfam_tophit_file):
            # verify checksum
            checksum_file = tigrfam_tophit_file + '.sha256'
            if os.path.exists(checksum_file):
              checksum = sha256(tigrfam_tophit_file)
              cur_checksum = open(checksum_file).readline().strip()
              if checksum == cur_checksum:
                continue

          tigrfam_file = os.path.join(assembly_dir, assembly_id + self.tigrfam_ext)
          if os.path.exists(tigrfam_file):
            tigrfam_files.append(tigrfam_file)

    print '  Number of unprocessed genomes: %d' % len(tigrfam_files)

    # populate worker queue with data to process
    workerQueue = mp.Queue()
    writerQueue = mp.Queue()

    for f in tigrfam_files:
      workerQueue.put(f)

    for _ in range(threads):
      workerQueue.put(None)

    try:
      workerProc = [mp.Process(target = self.__workerThread, args = (workerQueue, writerQueue)) for _ in range(threads)]
      writeProc = mp.Process(target = self.__writerThread, args = (len(tigrfam_files), writerQueue))

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
    p = TigrfamTopHit()
    p.run(args.genome_dir, args.threads)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
