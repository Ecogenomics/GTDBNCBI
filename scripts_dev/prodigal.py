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

__prog_name__ = 'prodigal.py'
__prog_desc__ = 'Run prodigal over a set of genomes.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
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
import multiprocessing as mp

from biolib.common import remove_extension
from biolib.external.execute import check_dependencies
from biolib.external.prodigal import Prodigal
from biolib.checksum import sha256


class RunProdigal(object):
    """Runs prodigal over a set of genomes.

    This script assumes that genomes are stored in individual
    directories for each user:

    <user_id>/<genome_id>/<genome_id>_genomic.fna
    <organism_dir>/<assembly_dir>/<genome_id>_genomic.fna

    where <genome_id>_genomic.fna is the genome file to be processed.
    """

    def __init__(self):
        check_dependencies(['prodigal'])

    def run(self, input_dir, tmp_dir, threads):
        # get path to all unprocessed genome files
        print 'Reading genomes.'
        genome_files = []
        for genome_dir in os.listdir(input_dir):
            cur_genome_dir = os.path.join(input_dir, genome_dir)
            if not os.path.isdir(cur_genome_dir):
                continue
              
            for assembly_id in os.listdir(cur_genome_dir):
                assembly_dir = os.path.join(cur_genome_dir, assembly_id)
                genome_id = assembly_id[0:assembly_id.find('_', 4)]

                # check if prodigal has already been called
                aa_gene_file = os.path.join(assembly_dir, 'prodigal', genome_id + '_protein.faa')
                if os.path.exists(aa_gene_file):
                    # verify checksum
                    checksum_file = aa_gene_file + '.sha256'
                    if os.path.exists(checksum_file):
                        checksum = sha256(aa_gene_file)
                        cur_checksum = open(checksum_file).readline().strip()
                        if checksum == cur_checksum:
                            continue

                genome_file = os.path.join(assembly_dir, assembly_id + '_genomic.fna')
                if os.path.exists(genome_file):
                    if os.stat(genome_file).st_size == 0:
                        print '[Warning] Genome file appears to be empty: %s' % genome_file
                    else:
                        genome_files.append(genome_file)

        print '  Number of unprocessed genomes: %d' % len(genome_files)

        # run prodigal on each genome
        print 'Running prodigal.'
        prodigal = Prodigal(cpus=threads)
        summary_stats = prodigal.run(genome_files, output_dir=tmp_dir)

        # move results into individual genome directories
        print 'Moving files and calculating checksums.'
        for genome_file in genome_files:
            genome_path, genome_id = ntpath.split(genome_file)
            genome_id = remove_extension(genome_id)
            
            aa_gene_file = os.path.join(tmp_dir, genome_id + '_genes.faa')
            nt_gene_file = os.path.join(tmp_dir, genome_id + '_genes.fna')
            gff_file = os.path.join(tmp_dir, genome_id + '.gff')

            genome_root = genome_id[0:genome_id.find('_', 4)]
            prodigal_path = os.path.join(genome_path, 'prodigal')
            if not os.path.exists(prodigal_path):
                os.makedirs(prodigal_path)
            new_aa_gene_file = os.path.join(prodigal_path, genome_root + '_protein.faa')
            new_nt_gene_file = os.path.join(prodigal_path, genome_root + '_protein.fna')
            new_gff_file = os.path.join(prodigal_path, genome_root + '_protein.gff')

            os.system('mv %s %s' % (aa_gene_file, new_aa_gene_file))
            os.system('mv %s %s' % (nt_gene_file, new_nt_gene_file))
            os.system('mv %s %s' % (gff_file, new_gff_file))

            # save translation table information
            translation_table_file = os.path.join(prodigal_path, 'prodigal_translation_table.tsv')
            fout = open(translation_table_file, 'w')
            fout.write('%s\t%d\n' % ('best_translation_table', summary_stats[genome_id].best_translation_table))
            fout.write('%s\t%.2f\n' % ('coding_density_4', summary_stats[genome_id].coding_density_4 * 100))
            fout.write('%s\t%.2f\n' % ('coding_density_11', summary_stats[genome_id].coding_density_11 * 100))
            fout.close()

            checksum = sha256(new_aa_gene_file)
            fout = open(new_aa_gene_file + '.sha256', 'w')
            fout.write(checksum)
            fout.close()

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('genome_dir', help='directory with genomes in individual directories')
  parser.add_argument('--tmp_dir', help='temporary directory for storing intermediate results', default='/tmp')
  parser.add_argument('-t', '--threads', type=int, help='number of threads', default=1)

  args = parser.parse_args()

  try:
    p = RunProdigal()
    p.run(args.genome_dir, args.tmp_dir, args.threads)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
