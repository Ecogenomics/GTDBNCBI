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

__prog_name__ = 'ncbi_genome_category.py'
__prog_desc__ = 'Identify genomes marked by NCBI as being a MAG or SAG.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2018'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import re
import argparse
from collections import defaultdict


class GenomeType(object):
    """Identify genomes marked by NCBI as being a MAG or SAG."""

    def __init__(self):
        pass

    def run(self, genbank_assembly_summary, refseq_assembly_summary, genome_file, output_file):
        """Identify genomes marked by NCBI as being a MAG or SAG."""
        
        # parse GenBank and RefSeq assembly summary report
        print('Parsing GenBank and RefSeq assembly reports.')
        gid_genome_type = {}
        for assembly_summary in [genbank_assembly_summary, refseq_assembly_summary]:
            with open(genbank_assembly_summary, encoding='utf-8') as f:
                for line in f:
                    if line.startswith('# assembly_accession'):
                        header = line.strip().split('\t')
                        exclude_idx = header.index('excluded_from_refseq')
                        continue
                    elif line.startswith('#'):
                        continue
                        
                    tokens = [t.strip() for t in line.split('\t')]
                    gid = tokens[0]
                    exclude = tokens[exclude_idx]
                    
                    if 'derived from single cell' in exclude:
                        gid_genome_type[gid] = 'SAG'
                    elif 'derived from metagenome' in exclude:
                        gid_genome_type[gid] = 'MAG'
                    elif 'derived from environmental source' in exclude:
                        gid_genome_type[gid] = 'ENV'
                    elif 'derived from surveillance project' in exclude:
                        pass
                    elif 'derived' in exclude:
                        # seems like an unhandled case
                        print('Unhandled case of derived genome: {}'.format(exclude))
                        sys.exit(-1)
        
        print(' - identified {:,} MAGs'.format(sum([1 for v in gid_genome_type.values() if v == 'MAG'])))
        print(' - identified {:,} SAGs'.format(sum([1 for v in gid_genome_type.values() if v == 'SAG'])))
        print(' - identified {:,} environmental genomes'.format(sum([1 for v in gid_genome_type.values() if v == 'ENV'])))

        # parse GBFF files
        print('Parsing GBFF file for each genome.')
        metagenome_pattern = re.compile(r'derived from(\w|\s)*metagenome', re.I)

        fout = open(output_file, 'w')
        fout.write('genome_id\tncbi_genome_category\tsource\n')
        fout_sanity_check = open(output_file + '.raw', 'w')
        genome_count = 0
        for idx, line in enumerate(open(genome_file)):
            if idx % 100 == 0:
                sys.stdout.write('==> Processed %d genomes.\r' % idx)
                sys.stdout.flush()

            gid, genome_dir, _fgid = line.strip().split('\t')

            types = set()

            if gid.startswith('U_'):
                source = 'user genome'
                types.add('MAG')
            elif gid in gid_genome_type:
                source = 'assembly report'
                types.add(gid_genome_type[gid])
            else:
                source = 'GBFF file'
                assembly_id = os.path.basename(os.path.normpath(genome_dir))
                wgs_file = os.path.join(genome_dir, assembly_id + '_genomic.gbff')
                
                for line in open(wgs_file):
                    if 'metagenome' in line:
                        if('/metagenome_source' in line
                                or re.search(metagenome_pattern, line)):
                            types.add('MAG')
                            fout_sanity_check.write('%s\t\%s' % (gid, line))
                    elif 'single cell' in line:
                        types.add('SAG')
                        fout_sanity_check.write('%s\t\%s' % (gid, line))
                    elif '/environmental_sample' in line or "derived from environmental source" in line:
                        types.add('ENV')
                        fout_sanity_check.write('%s\t\%s' % (gid, line))

            if 'MAG' in types and 'SAG' in types:
                print('[WARNING] Genome %s is annotated as both a MAG and SAG.' % gid)

            if 'MAG' in types:
                fout.write('%s\t%s\t%s\n' % (gid, 'derived from metagenome', source))
            elif 'SAG' in types:
                fout.write('%s\t%s\t%s\n' % (gid, 'derived from single cell', source))
            elif 'ENV' in types:
                fout.write('%s\t%s\t%s\n' %
                           (gid, 'derived from environmental sample', source))
        sys.stdout.flush()

        fout.close()
        fout_sanity_check.close()


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'genbank_assembly_summary', help='file from NCBI indicating metadata for genome assemblies in GenBank (ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt)')
    parser.add_argument(
        'refseq_assembly_summary', help='file from NCBI indicating metadata for genome assemblies in RefSeq (ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt)')
    parser.add_argument(
        'genome_file', help='file indicating path to GTDB genome files')
    parser.add_argument('output_file', help='output file')

    args = parser.parse_args()

    try:
        p = GenomeType()
        p.run(args.genbank_assembly_summary, 
                args.refseq_assembly_summary,
                args.genome_file, 
                args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
        raise
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
