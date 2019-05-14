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

__prog_name__ = 'rna_ltp.py'
__prog_desc__ = 'Identify, extract, and taxonomically classify 16S rRNA genes against the LTP DB.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2018'
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
import subprocess

from collections import defaultdict

from biolib.parallel import Parallel
from biolib.external.execute import check_dependencies

class RNA_LTP(object):
    """Identify, extract, and taxonomically classify 16S rRNA genes against the LTP DB."""

    def __init__(self):
        logger = logging.getLogger('')
        logger.setLevel(logging.DEBUG)
        log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                       datefmt="%Y-%m-%d %H:%M:%S")
        
        self.silva_output_dir = 'rna_silva_132' # needed to pick up previously identified and extracted 16S rRNA genes
        
        self.ltp_output_dir = 'rna_ltp_132'
        self.ltp_ssu_file = '/srv/whitlam/bio/db/silva/ltp/132/ltp_132.fna'
        self.ltp_taxonomy_file = '/srv/whitlam/bio/db/silva/ltp/132/ltp_132_taxonomy.tsv'

    def _producer(self, input_files):
        """Process each genome."""
        
        genome_file, ssu_file = input_files

        full_genome_dir, _ = ntpath.split(genome_file)

        output_dir = os.path.join(full_genome_dir, self.ltp_output_dir)

        # clean up old log files
        log_file = os.path.join(output_dir, 'genometk.log')
        if os.path.exists(log_file):
            os.remove(log_file)
          
        cmd_to_run = ['genometk', 
                        'rna', 
                        '--silent', 
                        '--cpus', 
                        '1', 
                        '--db', self.ltp_ssu_file, 
                        '--taxonomy_file', self.ltp_taxonomy_file,
                        '--rrna_file', ssu_file,
                        genome_file, 
                        'ssu', 
                        output_dir]
            
        proc = subprocess.Popen(cmd_to_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError("%r failed, status code %s stdout %r stderr %r" % (
                           cmd_to_run, proc.returncode, stdout, stderr))

        return output_dir

    def _progress(self, processed_items, total_items):
        return '  Processed %d of %d (%.2f%%) genomes.' % (processed_items,
                                                          total_items,
                                                          processed_items * 100.0 / total_items)

    def run(self, ncbi_genome_dir, user_genome_dir, cpus):
        """Create metadata by parsing assembly stats files."""

        input_files = []

        # generate metadata for NCBI assemblies
        if ncbi_genome_dir != 'NONE':
            print('Reading NCBI assembly directories.')
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
                        ssu_file = os.path.join(full_assembly_dir, self.silva_output_dir, 'ssu.fna')
                        if os.path.exists(ssu_file):
                            genome_file = os.path.join(full_assembly_dir, assembly_dir + '_genomic.fna')
                            input_files.append((genome_file, ssu_file))

        # generate metadata for user genomes
        if user_genome_dir != 'NONE':
            print('Reading user genome directories.')
            for user_id in os.listdir(user_genome_dir):
                full_user_dir = os.path.join(user_genome_dir, user_id)
                if not os.path.isdir(full_user_dir):
                    continue

                for genome_id in os.listdir(full_user_dir):
                    full_genome_dir = os.path.join(full_user_dir, genome_id)

                    ssu_file = os.path.join(full_assembly_dir, self.silva_output_dir, 'ssu.fna')
                    if os.path.exists(ssu_file):
                        genome_file = os.path.join(full_genome_dir, genome_id + '_genomic.fna')
                        input_files.append((genome_file, ssu_file))

                    print('Identified %d genomes to process.' % len(input_files))

        # process each genome
        print('Generating metadata for each genome:')
        parallel = Parallel(cpus = cpus)
        parallel.run(self._producer,
                      None,
                      input_files,
                      self._progress)

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ncbi_genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies or NONE to skip')
    parser.add_argument('user_genome_dir', help='base directory leading to user genomes or NONE to skip')
    parser.add_argument('-t', '--threads', help='number of CPUs to use', type=int, default=32)

    args = parser.parse_args()

    check_dependencies(['genometk'])

    try:
        p = RNA_LTP()
        p.run(args.ncbi_genome_dir, args.user_genome_dir, args.threads)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise