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

__prog_name__ = 'metadata_generate.py'
__prog_desc__ = 'Produce metadata files derived from nucleotide sequences and called genes.'

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


class Metadata(object):
  """Produce metadata files derived from nucleotide sequences and called genes."""

  def __init__(self):
    pass

  def _producer(self, input_files):
    """Process each genome."""

    genome_file, gff_file = input_files
    full_genome_dir, _ = ntpath.split(genome_file)

    # clean up old log files
    log_file = os.path.join(full_genome_dir, 'genometk.log')
    if os.path.exists(log_file):
      os.remove(log_file)

    # calculate metadata
    os.system('genometk nucleotide --silent %s %s' % (genome_file, full_genome_dir))
    os.system('genometk gene --silent %s %s %s' % (genome_file, gff_file, full_genome_dir))

    return full_genome_dir

  def _progress(self, processed_items, total_items):
      return '  Processed %d of %d (%.2f%%) genomes.' % (processed_items,
                                                          total_items,
                                                          processed_items * 100.0 / total_items)

  def run(self, ncbi_genome_dir, user_genome_dir, cpus):
    """Create metadata by parsing assembly stats files."""

    input_files = []

    # generate metadata for NCBI assemblies
    print 'Reading NCBI assembly directories.'
    processed_assemblies = defaultdict(list)
    for domain in ['archaea', 'bacteria']:
      domain_dir = os.path.join(ncbi_genome_dir, domain)
      for species_dir in os.listdir(domain_dir):
        full_species_dir = os.path.join(domain_dir, species_dir)
        for assembly_dir in os.listdir(full_species_dir):
          accession = assembly_dir[0:assembly_dir.find('_', 4)]

          processed_assemblies[accession].append(species_dir)
          if len(processed_assemblies[accession]) >= 2:
            continue

          full_assembly_dir = os.path.join(full_species_dir, assembly_dir)
          genome_file = os.path.join(full_assembly_dir, assembly_dir + '_genomic.fna')
          gff_file = os.path.join(full_assembly_dir, assembly_dir + '_genomic.gff')
          input_files.append([genome_file, gff_file])

    # generate metadata for user genomes
    print 'Reading user genome directories.'
    for user_id in os.listdir(user_genome_dir):
      full_user_dir = os.path.join(user_genome_dir, user_id)
      if not os.path.isdir(full_user_dir):
        continue

      for genome_id in os.listdir(full_user_dir):
        full_genome_dir = os.path.join(full_user_dir, genome_id)
        genome_file = os.path.join(full_genome_dir, genome_id + '_genomic.fna')
        gff_file = os.path.join(full_genome_dir, genome_id + '_protein.gff')
        input_files.append([genome_file, gff_file])

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
  parser.add_argument('user_genome_dir', help='base directory leading to user genomes')
  parser.add_argument('-c', '--cpus', help='number of CPUs to use', type=int, default=32)

  args = parser.parse_args()

  check_dependencies(['genometk'])

  try:
    p = Metadata()
    p.run(args.ncbi_genome_dir, args.user_genome_dir, args.cpus)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
