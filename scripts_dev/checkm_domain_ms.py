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

__prog_name__ = 'checkm_domain_ms.py'
__prog_desc__ = 'run CheckM domain-specific markers over a set of genomes'

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

  def run(self, genome_dir, domain, output_dir, cpus):
    """Applying CheckM to genomes."""

    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    tmp_dir = os.path.join(output_dir, 'genome_chunks')

    # determine gene files
    gene_files = []
    assembly_ids = {}
    for genome_id in os.listdir(genome_dir):
      cur_genome_dir = os.path.join(genome_dir, genome_id)
      if os.path.isdir(cur_genome_dir):
        for assembly_id in os.listdir(cur_genome_dir):
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

    # create simlinks to genomes in batches of 1000
    print 'Partitioning genomes into chunks of 1000.'
    num_chunks = 0
    for i, gene_file in enumerate(gene_files):
      if i % 1000 == 0:
        chunk_dir = os.path.join(tmp_dir, 'chunk%d' % num_chunks)
        os.makedirs(chunk_dir)
        num_chunks += 1

      os.system('ln -s %s %s' % (os.path.abspath(gene_file), os.path.join(chunk_dir, ntpath.basename(gene_file))))

    # apply CheckM to each set of 1000 genomes
    print 'Running CheckM on chunks:'
    for i in xrange(0, num_chunks):
      print '  Processing chunk %d of %d.' % (i+1, num_chunks)

      bin_dir = os.path.join(tmp_dir, 'chunk%d' % i)
      checkm_output_dir = os.path.join(output_dir, 'chunk%d' % i)
      if os.path.exists(checkm_output_dir):
        continue

      os.system('checkm taxonomy_wf --genes -x faa -t %d domain %s %s %s' % (cpus, domain, bin_dir, checkm_output_dir))

      qa_file = os.path.join(checkm_output_dir, 'qa.chunk%d.tsv' % i)
      os.system('checkm qa -t %d --tab_table -f %s %s %s' % (cpus, qa_file, os.path.join(checkm_output_dir, domain + '.ms'), checkm_output_dir))

    # create single file with CheckM results
    print 'Creating single file with CheckM results.'
    checkm_output = os.path.join(output_dir, 'checkm.profiles.domain_ms.tsv')
    fout = open(checkm_output, 'w')
    for i in xrange(0, num_chunks):
      profile_file = os.path.join(output_dir, 'chunk%d' % i, 'qa.chunk%d.tsv' % i)
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
  parser.add_argument('domain', help='domain level marker set to use', choices=['Bacteria', 'Archaea'])
  parser.add_argument('output_dir', help='output directory')
  parser.add_argument('-c', '--cpus', help='number of processors to use', type=int, default=16)

  args = parser.parse_args()

  try:
    runCheckm = RunCheckm()
    runCheckm.run(args.genome_dir, args.domain, args.output_dir, args.cpus)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
