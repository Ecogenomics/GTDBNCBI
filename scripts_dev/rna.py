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

__prog_name__ = 'rna.py'
__prog_desc__ = 'Identify, extract, and taxonomically classify rRNA genes in genomes.'

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
import logging
import ntpath
import argparse
from collections import defaultdict

from biolib.parallel import Parallel
from biolib.external.execute import check_dependencies

class RNA(object):
  """Identify, extract, and taxonomically classify rRNA genes in genomes."""

  def __init__(self):
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                   datefmt="%Y-%m-%d %H:%M:%S")

  def _producer(self, genome_file):
    """Process each genome."""

    full_genome_dir, _ = ntpath.split(genome_file)

    output_dir = os.path.join(full_genome_dir, self.output_dir)

    # clean up old log files
    log_file = os.path.join(output_dir, 'genometk.log')
    if os.path.exists(log_file):
      os.remove(log_file)

    os.system('genometk rna --silent --cpus 1 --db %s --taxonomy_file %s %s %s %s' % (self.db, self.taxonomy, genome_file, self.rna_gene, output_dir))

    return output_dir

  def _progress(self, processed_items, total_items):
      return '  Processed %d of %d (%.2f%%) genomes.' % (processed_items,
                                                          total_items,
                                                          processed_items * 100.0 / total_items)

  def run(self, rna_gene, ncbi_genome_dir, user_genome_dir, cpus):
    """Create metadata by parsing assembly stats files."""
    
    print 'Running with GreenGenes database'
    if rna_gene == 'ssu':
       self._run(rna_gene, ncbi_genome_dir, user_genome_dir, 'GG', cpus)
    
    print 'Running with SILVA database'
    self._run(rna_gene, ncbi_genome_dir, user_genome_dir, 'SILVA', cpus)
    
  def _run(self, rna_gene, ncbi_genome_dir, user_genome_dir, ssu_db, cpus):
    """Create metadata by parsing assembly stats files."""
    
    self.rna_gene = rna_gene
    
    if ssu_db == 'GG':
        # Greengenes data files and desired output
        if rna_gene == 'ssu':
            self.db = '/srv/db/gg/2013_08/gg_13_8_otus/rep_set/99_otus.fasta'
            self.taxonomy = '/srv/db/gg/2013_08/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt'
            self.output_dir = 'ssu_gg'
        else:
            print 'There is no LSU database for GG.'
            sys.exit()
    elif ssu_db == 'SILVA':
        # Silva info
        if rna_gene == 'ssu':
            self.db = '/srv/whitlam/bio/db/silva/123.1/SILVA_123.1_SSURef_Nr99_tax_silva.fasta'
            self.taxonomy = '/srv/whitlam/bio/db/silva/123.1/silva_taxonomy.ssu.tsv'
        elif rna_gene == 'lsu_23S':
            self.db = '/srv/db/silva/123.1/SILVA_123.1_LSURef_tax_silva.fasta'
            self.taxonomy = '/srv/whitlam/bio/db/silva/123.1/silva_taxonomy.lsu.tsv'
        self.output_dir = 'rna_silva'

    input_files = []

    # generate metadata for NCBI assemblies
    print 'Reading NCBI assembly directories.'
    processed_assemblies = defaultdict(list)
    for domain in ['archaea', 'bacteria']:
      domain_dir = os.path.join(ncbi_genome_dir, domain)
      if not os.path.exists(domain_dir):
        continue

      for species_dir in os.listdir(domain_dir):
        full_species_dir = os.path.join(domain_dir, species_dir)
        for assembly_dir in os.listdir(full_species_dir):
          accession = assembly_dir[0:assembly_dir.find('_', 4)]

          processed_assemblies[accession].append(species_dir)
          if len(processed_assemblies[accession]) >= 2:
            continue

          full_assembly_dir = os.path.join(full_species_dir, assembly_dir)

          #***if os.path.exists(os.path.join(full_assembly_dir, self.output_dir)):
          #***  continue

          genome_file = os.path.join(full_assembly_dir, assembly_dir + '_genomic.fna')
          input_files.append(genome_file)

    # generate metadata for user genomes
    if user_genome_dir != 'NONE':
        print 'Reading user genome directories.'
        for user_id in os.listdir(user_genome_dir):
          full_user_dir = os.path.join(user_genome_dir, user_id)
          if not os.path.isdir(full_user_dir):
            continue

          for genome_id in os.listdir(full_user_dir):
            full_genome_dir = os.path.join(full_user_dir, genome_id)

            if os.path.exists(os.path.join(full_genome_dir, self.output_dir)):
                continue

            genome_file = os.path.join(full_genome_dir, genome_id + '_genomic.fna')
            input_files.append(genome_file)

        print 'Identified %d genomes to process.' % len(input_files)

    # process each genome
    print 'Generating metadata for each genome:'
    parallel = Parallel(cpus = cpus)
    parallel.run(self._producer,
                  None,
                  input_files,
                  self._progress)

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('rna_gene', choices=['ssu', 'lsu_23S'], help="rRNA gene to process")
  parser.add_argument('ncbi_genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies')
  parser.add_argument('user_genome_dir', help='base directory leading to user genomes or NONE to skip')
  parser.add_argument('-t', '--threads', help='number of CPUs to use', type=int, default=32)

  args = parser.parse_args()

  check_dependencies(['genometk'])

  try:
    p = RNA()
    p.run(args.rna_gene, args.ncbi_genome_dir, args.user_genome_dir, args.threads)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
