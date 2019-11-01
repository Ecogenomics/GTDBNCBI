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

__prog_name__ = 'ncbi_assembly_file_metadata.py'
__prog_desc__ = 'Produce filtered metadata file from NCBI assembly metadata file.'

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
import argparse
from collections import defaultdict

from numpy import (zeros as np_zeros)


class Metadata(object):
    """Create metadata file from the assembly stats file of each NCBI assembly."""

    def __init__(self):
        self.fields = {'bioproject': 'ncbi_bioproject',
                       'wgs_master': 'ncbi_wgs_master',
                       'refseq_category': 'ncbi_refseq_category',
                       'species_taxid': 'ncbi_species_taxid',
                       'isolate': 'ncbi_isolate',
                       'version_status': 'ncbi_version_status',
                       'seq_rel_date': 'ncbi_seq_rel_date',
                       'asm_name': 'ncbi_asm_name',
                       'gbrs_paired_asm': 'ncbi_gbrs_paired_asm',
                       'paired_asm_comp': 'ncbi_paired_asm_comp',
                       'relation_to_type_material': 'ncbi_type_material_designation'}

    def run(self, refseq_bacteria_assembly_summary_file,
            refseq_archaea_assembly_summary_file,
            genbank_bacteria_assembly_summary_file,
            genbank_archaea_assembly_summary_file, genome_id_file, output_file):
        """Create metadata by parsing NCBI assembly metadata file."""

        # get identifier of genomes in GTDB
        genome_ids = set()
        for line in open(genome_id_file):
            if line[0] == '#':
                continue

            if '\t' in line:
                genome_id = line.strip().split('\t')[0]
            else:
                genome_id = line.strip().split(',')[0]

            if genome_id.startswith('GCA_'):
                genome_id = 'GB_' + genome_id
            elif genome_id.startswith('GCF_'):
                genome_id = 'RS_' + genome_id

            genome_ids.add(genome_id)

        # write out metadata
        fout = open(output_file, 'w')
        fout.write('genome_id')

        write_header = True
        indice_wgs = None
        for assembly_file in [refseq_bacteria_assembly_summary_file,
                              refseq_archaea_assembly_summary_file,
                              genbank_bacteria_assembly_summary_file,
                              genbank_archaea_assembly_summary_file]:
            with open(assembly_file) as f:
                f.readline()  # first comment line
                headers = f.readline().rstrip().split('\t')

                indices = []
                for i, header in enumerate(headers):
                    if header in self.fields:
                        if write_header:
                            fout.write('\t' + self.fields[header])
                        indices.append(i)
                        if write_header and header == 'wgs_master':
                            fout.write('\t' + 'ncbi_wgs_formatted')
                            indice_wgs = i

                if write_header:
                    fout.write('\n')
                    write_header = False

                print indice_wgs

                for line in f:
                    line_split = line.rstrip('\n').split('\t')

                    genome_id = line_split[0]
                    if genome_id.startswith('GCA_'):
                        genome_id = 'GB_' + genome_id
                    elif genome_id.startswith('GCF_'):
                        genome_id = 'RS_' + genome_id

                    if genome_id in genome_ids:
                        fout.write(genome_id)
                        for i in indices:
                            if indice_wgs == i:
                                fout.write('\t' + line_split[i])
                                fout.write(
                                    '\t' + self.format_wgs(line_split[i]))
                            else:
                                fout.write('\t' + line_split[i])
                        fout.write('\n')

        fout.close()

    def format_wgs(self, wgs_accession):
        if not wgs_accession:
            return ""
        wgs_acc, version = wgs_accession.split('.')
        idx = [ch.isdigit() for ch in wgs_acc].index(True)
        wgs_id = wgs_acc[0:idx] + str(version).zfill(2)
        return wgs_id


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('refseq_bacteria_assembly_summary_file',
                        help='RefSeq assembly summary file (assembly_summary_refseq.txt)')
    parser.add_argument('refseq_archaea_assembly_summary_file',
                        help='RefSeq assembly summary file (assembly_summary_refseq.txt)')
    parser.add_argument('genbank_bacteria_assembly_summary_file',
                        help='GenBank assembly summary file (assembly_summary_genbank.txt)')
    parser.add_argument('genbank_archaea_assembly_summary_file',
                        help='GenBank assembly summary file (assembly_summary_genbank.txt)')
    parser.add_argument(
        'genome_id_file', help='genome identifiers for genomes in GTDB')
    parser.add_argument('output_file', help='output metadata file')

    args = parser.parse_args()

    try:
        p = Metadata()
        p.run(args.refseq_bacteria_assembly_summary_file,
              args.refseq_archaea_assembly_summary_file,
              args.genbank_bacteria_assembly_summary_file,
              args.genbank_archaea_assembly_summary_file, args.genome_id_file, args.output_file)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
