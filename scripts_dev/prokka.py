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

__prog_name__ = 'prokka.py'
__prog_desc__ = 'Run Prokka over a set of genomes.'

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
from shutil import copyfile
import multiprocessing as mp

from biolib.checksum import sha256
from biolib.external.execute import check_dependencies


class Prokka(object):
  """Runs Prokka over a set of genomes.

  This script assumes that genomes are stored in individual
  directories in the following format:

  <domain>/<genome_id>/<assembly_id>/<assembly_id>_protein.faa

  where <domain> is either 'archaea' or 'bacteria', and there
    may be multiple assembly_id for a given genome_id. These
    typically represent different strains from a species.

  This is the directory structure which results from extract_ncbi.py.
  """

  def __init__(self):
    check_dependencies(['prokka'])

    self.genome_file_ext = '_genomic.fna'
    self.protein_aa_file_ext = '_protein.faa'
    self.protein_nt_file_ext = '_protein.fna'

  def __workerThread(self, domain, queueIn, queueOut):
    """Process each data item in parallel."""
    while True:
      genome_file = queueIn.get(block=True, timeout=None)
      if genome_file == None:
        break

      assembly_dir, filename = os.path.split(genome_file)
      prefix = filename.replace(self.genome_file_ext, '')
      output_dir = os.path.join(assembly_dir, 'prokka')
      if os.path.exists(output_dir):
        queueOut.put(genome_file)
        continue
        
      os.makedirs(output_dir)
      prokka_out = os.path.join(output_dir, 'prokka.out')

      cmd = 'prokka --force --kingdom %s --prefix %s --outdir %s --cpus 1 %s 2> %s' % (domain, prefix, output_dir, genome_file, prokka_out)
      os.system(cmd)

      # calculate checksum
      prokka_gene_file = os.path.join(output_dir, prefix + '.faa')
      checksum = sha256(prokka_gene_file)
      fout = open(prokka_gene_file + '.sha256', 'w')
      fout.write(checksum)
      fout.close()
      
      # copy files
      protein_nt_file = os.path.join(output_dir, prefix + '.ffn')
      output_file = os.path.join(assembly_dir, prefix + self.protein_nt_file_ext)
      copyfile(protein_nt_file, output_file)
      
      protein_aa_file = os.path.join(output_dir, prefix + '.faa')
      output_file = os.path.join(assembly_dir, prefix + self.protein_aa_file_ext)
      copyfile(protein_aa_file, output_file)
      
      # allow results to be processed or written to file
      queueOut.put(genome_file)

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

  def run(self, genome_dir, domain, threads):
    # get path to all unprocessed genome gene files
    print 'Reading genomes.'
    genome_files = []
    for genome_id in os.listdir(genome_dir):
      cur_genome_dir = os.path.join(genome_dir, genome_id)
      if os.path.isdir(cur_genome_dir):
        for assembly_id in os.listdir(cur_genome_dir):
          assembly_dir = os.path.join(cur_genome_dir, assembly_id)
          
          protein_file = os.path.join(assembly_dir, assembly_id + self.protein_aa_file_ext)
          if os.path.exists(protein_file):
            continue

          prokka_dir = os.path.join(assembly_dir, 'prokka')
          if os.path.exists(prokka_dir):
            continue
            
          prokka_file = os.path.join(prokka_dir, assembly_id + '.faa')
          if os.path.exists(prokka_file):
            # verify checksum
            checksum_file = prokka_file + '.sha256'
            if os.path.exists(checksum_file):
              checksum = sha256(prokka_file)
              cur_checksum = open(checksum_file).readline().strip()
              if checksum == cur_checksum:
                continue

          genome_file = os.path.join(assembly_dir, assembly_id + self.genome_file_ext)
          if os.path.exists(genome_file):
            genome_files.append(genome_file)

    print '  Number of unprocessed genomes: %d\n' % len(genome_files)

    # populate worker queue with data to process
    workerQueue = mp.Queue()
    writerQueue = mp.Queue()

    for f in genome_files:
      workerQueue.put(f)

    for _ in range(threads):
      workerQueue.put(None)

    try:
      workerProc = [mp.Process(target = self.__workerThread, args = (domain, workerQueue, writerQueue)) for _ in range(threads)]
      writeProc = mp.Process(target = self.__writerThread, args = (len(genome_files), writerQueue))

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
  parser.add_argument('domain', help='domain to use for annotation by Prokka (e.g., Archaea, Bacteria)')
  parser.add_argument('-t', '--threads', type=int, help='number of threads', default=1)

  args = parser.parse_args()

  try:
    p = Prokka()
    p.run(args.genome_dir, args.domain, args.threads)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
