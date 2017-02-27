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

__prog_name__ = 'tigrfam.py'
__prog_desc__ = 'Run the TIGRfam HMMs over a set of genomes.'

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


class Tigrfam(object):
  """Runs TIGRfam HMMs over a set of genomes.

  This script assumes that genomes are stored in individual
  directories in the following format:

  <domain>/<genome_id>/<assembly_id>/<assembly_id>_protein.faa

  where <domain> is either 'archaea' or 'bacteria', and there
    may be multiple assembly_id for a given genome_id. These
    typically represent different strains from a species.

  This is the directory structure which results from extract_ncbi.py.
  """

  def __init__(self):
    check_dependencies(['hmmsearch'])

    self.tigrfam_hmms = '/srv/whitlam/bio/db/tigrfam/15.0/TIGRFAMs_15.0_HMM/tigrfam.hmm'

    self.protein_file_ext = '_protein.faa'

  def __workerThread(self, queueIn, queueOut):
    """Process each data item in parallel."""
    while True:
      gene_file = queueIn.get(block=True, timeout=None)
      if gene_file == None:
        break

      assembly_dir, filename = os.path.split(gene_file)
      
      running_file = os.path.join(assembly_dir, filename.replace(self.protein_file_ext, '_tigrfam.running'))
      if not os.path.exists(running_file):
        fout = open(running_file, 'w')
        fout.write('running')
        fout.close()
        
        output_hit_file = os.path.join(assembly_dir, filename.replace(self.protein_file_ext, '_tigrfam.tsv'))
        hmmsearch_out = os.path.join(assembly_dir, filename.replace(self.protein_file_ext, '_tigrfam.out'))
        cmd = 'hmmsearch -o %s --tblout %s --noali --notextw --cut_nc --cpu 1 %s %s' % (hmmsearch_out, output_hit_file, self.tigrfam_hmms, gene_file)
        os.system(cmd)

        # calculate checksum
        checksum = sha256(output_hit_file)
        fout = open(output_hit_file + '.sha256', 'w')
        fout.write(checksum)
        fout.close()

        if os.path.exists(running_file):
            os.remove(running_file)

      # allow results to be processed or written to file
      queueOut.put(gene_file)

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

  def run(self, genome_dir, genome_report, threads):
    # get list of genomes to consider
    genomes_to_consider = set()
    for line in open(genome_report):
        line_split = line.strip().split('\t')
        genome_id = line_split[1]
        
        attributes = line_split[2].split(';')
        for attribute in attributes:
            if attribute == 'new' or attribute == 'modified':
                genomes_to_consider.add(genome_id)

    print 'Identified %d genomes as new or modified.' % len(genomes_to_consider)
    
    # get path to all unprocessed genome gene files
    print 'Reading genomes.'
    gene_files = []
    for species_dir in os.listdir(genome_dir):
      cur_genome_dir = os.path.join(genome_dir, species_dir)
      if os.path.isdir(cur_genome_dir):
        for assembly_id in os.listdir(cur_genome_dir):
          prodigal_dir = os.path.join(cur_genome_dir, assembly_id, 'prodigal')
          genome_id = assembly_id[0:assembly_id.find('_', 4)]

          tigrfam_file = os.path.join(prodigal_dir, genome_id + '_tigrfam.tsv')
          if os.path.exists(tigrfam_file):
            # verify checksum
            checksum_file = tigrfam_file + '.sha256'
            if os.path.exists(checksum_file):
              checksum = sha256(tigrfam_file)
              cur_checksum = open(checksum_file).readline().strip()
              if checksum == cur_checksum:
                if genome_id in genomes_to_consider:
                    print '[WARNING] Genome %s is marked as new or modified, but already has TIGRfam annotations.' % genome_id
                    print '[WARNING] Genome is being skipped!'
                continue
                
            print '[WARNING] Genome %s has TIGRfam annotations, but an invalid checksum and was not marked for reannotation.' % genome_id
            print '[WARNING] Genome will be reannotated.'
      
          elif genome_id not in genomes_to_consider:
            print '[WARNING] Genome %s has no TIGRfam annotations, but is also not marked for processing?' % genome_id
            print '[WARNING] Genome will be reannotated!'

          gene_file = os.path.join(prodigal_dir, genome_id + self.protein_file_ext)
          if os.path.exists(gene_file):
            if os.stat(gene_file).st_size == 0:
                print '[Warning] Protein file appears to be empty: %s' % gene_file
            else:
                gene_files.append(gene_file)

    print '  Number of unprocessed genomes: %d' % len(gene_files)

    # populate worker queue with data to process
    workerQueue = mp.Queue()
    writerQueue = mp.Queue()

    for f in gene_files:
      workerQueue.put(f)

    for _ in range(threads):
      workerQueue.put(None)

    try:
      workerProc = [mp.Process(target = self.__workerThread, args = (workerQueue, writerQueue)) for _ in range(threads)]
      writeProc = mp.Process(target = self.__writerThread, args = (len(gene_files), writerQueue))

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
  parser.add_argument('genome_report', help='report log indicating new, modified, unmodified, ..., genomes')
  parser.add_argument('-t', '--threads', type=int, help='number of threads', default=1)

  args = parser.parse_args()

  try:
    p = Tigrfam()
    p.run(args.genome_dir, args.genome_report, args.threads)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
