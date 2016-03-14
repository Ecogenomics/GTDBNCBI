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

__prog_name__ = 'checkm_lineage_ms.py'
__prog_desc__ = 'run CheckM lineage-specific markers over a set of genomes'

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
import tempfile
import ntpath
from string import maketrans

import biolib.seq_io as seq_io

class RunCheckm(object):
  """Apply CheckM to a large set of genomes.

  This script assumes that genomes are stored in individual
  directories in the following format:

  <domain>/<genome_id>/<assembly_id>/<assembly_id>_protein.faa

  where <domain> is either 'archaea' or 'bacteria', and there
    may be multiple assembly_id for a given genome_id. These
    typically represent different strains from a species.

  This is the directory structure which results from extract_ncbi.py.
  """

  def __init__(self):
    """Initialization."""
    pass

  def run(self, genome_dir, genome_report, cpus, output_dir):
    """Applying CheckM to genomes."""

    tmp_dir = os.path.join(output_dir, 'genome_chunks')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
        # get list of genomes to consider
        genomes_to_consider = set()
        for line in open(genome_report):
            line_split = line.strip().split('\t')
            if line_split[2] == 'new' or line_split[2] == 'modified':
                genomes_to_consider.add(line_split[1])
                
        print 'Identified %d genomes as new or modified.' % len(genomes_to_consider)
                
        # determine gene files
        gene_files = []
        assembly_ids = {}
        for genome_id in os.listdir(genome_dir):
          cur_genome_dir = os.path.join(genome_dir, genome_id)
          if os.path.isdir(cur_genome_dir):
            for assembly_id in os.listdir(cur_genome_dir):
              genome_id = assembly_id[0:assembly_id.find('_', 4)]
              if genome_id not in genomes_to_consider:
                continue
                
              assembly_dir = os.path.join(cur_genome_dir, assembly_id)

              if assembly_id in assembly_ids:
                print 'Duplicate assembly_id:'
                print assembly_ids[assembly_id]
                print assembly_dir
                continue

              assembly_ids[assembly_id] = assembly_dir

              gene_file = os.path.join(assembly_dir, assembly_id + '_protein.faa')
              if os.path.exists(gene_file):
                gene_files.append(gene_file)

        print '  Identified %d gene files.' % len(gene_files)

        # copy genomes in batches of 1000
        ambiguous_aa = maketrans('BJZ', 'XXX')
        
        print 'Partitioning genomes into chunks of 1000.'
        num_chunks = 0
        for i, gene_file in enumerate(gene_files):
          if i % 1000 == 0:
            chunk_dir = os.path.join(tmp_dir, 'chunk%d' % num_chunks)
            print chunk_dir
            os.makedirs(chunk_dir)
            num_chunks += 1

          # copy sequences and replace ambiguous bases with an X (unknown)
          # as pplacer is not compatible with these characters
          fout = open(os.path.join(chunk_dir, ntpath.basename(gene_file)), 'w')
          for seq_id, seq, annotation in  seq_io.read_seq(gene_file, True):
              fout.write('>' + seq_id + ' ' + annotation + '\n')
              fout.write(seq.translate(ambiguous_aa) + '\n')
    else:
        # just determine number of "chunk" directories
        num_chunks = 0
        for d in os.listdir(tmp_dir):
            if 'chunk' in d:
                num_chunks += 1

    # apply CheckM to each set of 1000 genomes
    print 'Running CheckM on chunks:'
    for i in xrange(0, num_chunks):
      print '  Processing chunk %d of %d.' % (i+1, num_chunks)

      bin_dir = os.path.join(tmp_dir, 'chunk%d' % i)
      checkm_output_dir = os.path.join(output_dir, 'chunk%d' % i)
      if os.path.exists(checkm_output_dir):
        continue

      os.makedirs(checkm_output_dir)
      os.system('checkm lineage_wf --pplacer_threads 10 --genes -x faa -t %d %s %s' % (cpus, bin_dir, checkm_output_dir))

      tree_qa_file = os.path.join(checkm_output_dir, 'tree_qa.o2.chunk%d.tsv' % i)
      os.system('checkm tree_qa -o 2 --tab_table -f %s %s' % (tree_qa_file, checkm_output_dir))

      qa_file = os.path.join(checkm_output_dir, 'qa.chunk%d.tsv' % i)
      os.system('checkm qa -t %d --tab_table -f %s %s %s' % (cpus, qa_file, os.path.join(checkm_output_dir, 'lineage.ms'), checkm_output_dir))

      profile_file = os.path.join(checkm_output_dir, 'profile.chunk%d.tsv' % i)
      os.system('checkm join_tables -f %s %s %s' % (profile_file, qa_file, tree_qa_file))

    # create single file with CheckM results
    print 'Creating single file with CheckM results.'
    checkm_output = os.path.join(output_dir, 'checkm.profiles.tsv')
    fout = open(checkm_output, 'w')
    for i in xrange(0, num_chunks):
      profile_file = os.path.join(output_dir, 'chunk%d' % i, 'profile.chunk%d.tsv' % i)
      with open(profile_file) as f:
        if i != 0:
          f.readline()

        for line in f:
          fout.write(line)
    fout.close()

    print '  CheckM results written to: %s' % checkm_output

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('genome_dir', help='directory containing genomes in individual directories')
  parser.add_argument('genome_report', help='report log indicating new, modified, unmodified, ..., genomes')
  parser.add_argument('output_dir', help='output directory')
  parser.add_argument('-c', '--cpus', help='number of processors to use', type=int, default=16)

  args = parser.parse_args()

  try:
    runCheckm = RunCheckm()
    runCheckm.run(args.genome_dir, args.genome_report, args.cpus, args.output_dir)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
