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

__prog_name__ = 'rna_silva.py'
__prog_desc__ = 'Identify, extract, and taxonomically classify rRNA genes in genomes against SILVA DB.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.3'
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


class RNA(object):
    """Identify, extract, and taxonomically classify rRNA genes in genomes."""

    def __init__(self):
        logger = logging.getLogger('')
        logger.setLevel(logging.DEBUG)
        log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                       datefmt="%Y-%m-%d %H:%M:%S")

        self.silva_output_dir = 'rna_silva_132'

        self.silva_ssu_ref_file = '/srv/whitlam/bio/db/silva/132/SILVA_132_SSURef_Nr99_tax_silva.fasta'
        self.silva_lsu_ref_file = '/srv/whitlam/bio/db/silva/132/SILVA_132_LSURef_tax_silva.fasta'

        self.silva_ssu_taxonomy_file = '/srv/whitlam/bio/db/silva/132/silva_taxonomy.ssu.tsv'
        self.silva_lsu_taxonomy_file = '/srv/whitlam/bio/db/silva/132/silva_taxonomy.lsu.tsv'

    def _producer(self, genome_file):
        """Process each genome."""

        full_genome_dir, _ = ntpath.split(genome_file)

        output_dir = os.path.join(full_genome_dir, self.output_dir)

        # clean up old log files
        log_file = os.path.join(output_dir, 'genometk.log')
        if os.path.exists(log_file):
            os.remove(log_file)

        if self.rna_gene == 'lsu_5S':
            cmd_to_run = ['genometk', 'rna', '--silent', '--cpus', '1',
                          '--min_len', '80', genome_file, self.rna_gene, output_dir]
        else:
            cmd_to_run = ['genometk', 'rna', '--silent', '--cpus', '1', '--db', self.db,
                          '--taxonomy_file', self.taxonomy, genome_file, self.rna_gene, output_dir]

        proc = subprocess.Popen(
            cmd_to_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError("%r failed, status code %s stdout %r stderr %r" % (
                cmd_to_run, proc.returncode, stdout, stderr))

        return output_dir

    def _progress(self, processed_items, total_items):
        return '  Processed %d of %d (%.2f%%) genomes.' % (processed_items,
                                                           total_items,
                                                           processed_items * 100.0 / total_items)

    def run(self, rna_gene, ncbi_genome_dir, user_genome_dir, cpus):
        """Create metadata by parsing assembly stats files."""

        # GreenGenes is effectively dead at this point
        #***print 'Running with GreenGenes database'
        #***if rna_gene == 'ssu':
        #***   self._run(rna_gene, ncbi_genome_dir, user_genome_dir, 'GG', cpus)

        print('Running with SILVA database')
        self._run(rna_gene, ncbi_genome_dir, user_genome_dir, 'SILVA', cpus)

    def _run(self, rna_gene, ncbi_genome_dir, user_genome_dir, db, cpus):
        """Create metadata by parsing assembly stats files."""

        self.rna_gene = rna_gene

        if db == 'SILVA':
            # Silva info
            if rna_gene == 'ssu':
                self.db = self.silva_ssu_ref_file
                self.taxonomy = self.silva_ssu_taxonomy_file
            elif rna_gene == 'lsu_23S':
                self.db = self.silva_lsu_ref_file
                self.taxonomy = self.silva_lsu_taxonomy_file
            elif rna_gene == 'lsu_5S':
                print('We currently do not curate against a 5S database, but do identify these sequences for quality assessment purposes.')
            self.output_dir = self.silva_output_dir
        else:
            print('Unrecognized database: %s' % db)
            sys.exit(-1)

        input_files = []

        # generate metadata for NCBI assemblies
        if ncbi_genome_dir != 'NONE':
            print('Reading NCBI assembly directories.')
            processed_assemblies = defaultdict(list)
            rfq_dir = os.path.join(ncbi_genome_dir, 'refseq', 'GCF')
            gbk_dir = os.path.join(ncbi_genome_dir, 'genbank', 'GCA')

            for input_dir in (rfq_dir, gbk_dir):
                for first_three in os.listdir(input_dir):
                    onethird_species_dir = os.path.join(input_dir, first_three)
                    print(onethird_species_dir)
                    if os.path.isfile(onethird_species_dir):
                        continue
                    for second_three in os.listdir(onethird_species_dir):
                        twothird_species_dir = os.path.join(
                            onethird_species_dir, second_three)
                        # print twothird_species_dir
                        if os.path.isfile(twothird_species_dir):
                            continue
                        for third_three in os.listdir(twothird_species_dir):
                            threethird_species_dir = os.path.join(
                                twothird_species_dir, third_three)
                            # print threethird_species_dir
                            if os.path.isfile(threethird_species_dir):
                                continue
                            for complete_name in os.listdir(threethird_species_dir):
                                assembly_dir = os.path.join(
                                    threethird_species_dir, complete_name)
                                if os.path.isfile(assembly_dir):
                                    continue

                                accession = complete_name[0:complete_name.find(
                                    '_', 4)]

                                processed_assemblies[accession].append(
                                    assembly_dir)
                                if len(processed_assemblies[accession]) >= 2:
                                    continue

                                hmm_results_file = os.path.join(
                                    assembly_dir, self.output_dir, rna_gene + '.hmm_summary.tsv')
                                if os.path.exists(hmm_results_file):
                                    continue

                                genome_file = os.path.join(
                                    assembly_dir, complete_name + '_genomic.fna')
                                input_files.append(genome_file)

        # generate metadata for user genomes
        if user_genome_dir != 'NONE':
            print('Reading user genome directories.')
            for user_id in os.listdir(user_genome_dir):
                full_user_dir = os.path.join(user_genome_dir, user_id)
                if not os.path.isdir(full_user_dir):
                    continue

                for genome_id in os.listdir(full_user_dir):
                    full_genome_dir = os.path.join(full_user_dir, genome_id)

                    hmm_results_file = os.path.join(
                        full_genome_dir, self.output_dir, rna_gene + '.hmm_summary.tsv')
                    if os.path.exists(hmm_results_file):
                        continue

                    genome_file = os.path.join(
                        full_genome_dir, genome_id + '_genomic.fna')
                    input_files.append(genome_file)

            print('Identified %d genomes to process.' % len(input_files))

        # process each genome
        print('Generating metadata for each genome:')
        parallel = Parallel(cpus=cpus)
        parallel.run(self._producer,
                     None,
                     input_files,
                     self._progress)


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('rna_gene', choices=[
                        'ssu', 'lsu_23S', 'lsu_5S'], help="rRNA gene to process")
    parser.add_argument(
        'ncbi_genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies or NONE to skip')
    parser.add_argument(
        'user_genome_dir', help='base directory leading to user genomes or NONE to skip')
    parser.add_argument('-t', '--threads',
                        help='number of CPUs to use', type=int, default=32)

    args = parser.parse_args()

    check_dependencies(['genometk'])

    try:
        p = RNA()
        p.run(args.rna_gene, args.ncbi_genome_dir,
              args.user_genome_dir, args.threads)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
