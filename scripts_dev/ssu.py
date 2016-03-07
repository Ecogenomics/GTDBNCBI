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

__prog_name__ = 'ssu.py'
__prog_desc__ = 'Identify, extract, and taxonomically classify 16S rRNA genes in genomes.'

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
import ntpath
import argparse
from collections import defaultdict

from biolib.parallel import Parallel
from biolib.external.execute import check_dependencies

class SSU(object):
  """Identify, extract, and taxonomically classify 16S rRNA genes in genomes."""

  def __init__(self):
    pass

  def _producer(self, genome_file):
    """Process each genome."""

    full_genome_dir, _ = ntpath.split(genome_file)

    output_dir = os.path.join(full_genome_dir, self.output_dir)

    # clean up old log files
    log_file = os.path.join(output_dir, 'genometk.log')
    if os.path.exists(log_file):
      os.remove(log_file)

    os.system('genometk ssu --silent --cpus 1 %s %s %s %s' % (genome_file, self.ssu_db, self.ssu_taxonomy, output_dir))

    return output_dir

  def _progress(self, processed_items, total_items):
      return '  Processed %d of %d (%.2f%%) genomes.' % (processed_items,
                                                          total_items,
                                                          processed_items * 100.0 / total_items)

  def run(self, ncbi_genome_dir, user_genome_dir, ssu_db, cpus):
    """Create metadata by parsing assembly stats files."""
    
    if ssu_db == 'GG':
        # Greengenes data files and desired output
        self.ssu_db = '/srv/db/gg/2013_08/gg_13_8_otus/rep_set/99_otus.fasta'
        self.ssu_taxonomy = '/srv/db/gg/2013_08/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt'
        self.output_dir = 'ssu_gg_2013_08'
    elif ssu_db == 'SILVA':
        # Silva info
        self.ssu_db = '/srv/whitlam/bio/db/silva/119/SILVA_119.pintail_100.ncbi_ids.fna'
        self.ssu_taxonomy = '/srv/whitlam/bio/db/silva/119/ssuref_119_maria.2015-09-28.filled.ncbi_ids.tsv'
        self.output_dir = 'ssu_silva_199_gg_taxonomy'

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

          if os.path.exists(os.path.join(full_assembly_dir, self.output_dir)):
            continue

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
  parser.add_argument('ncbi_genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies')
  parser.add_argument('user_genome_dir', help='base directory leading to user genomes or NONE to skip')
  parser.add_argument('ssu_db', choices=['GG', 'SILVA'], help='SSU database to use for assigning taxonomy')
  
  parser.add_argument('-t', '--threads', help='number of CPUs to use', type=int, default=32)

  args = parser.parse_args()

  check_dependencies(['genometk'])

  try:
    p = SSU()
    p.run(args.ncbi_genome_dir, args.user_genome_dir, args.ssu_db, args.threads)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
