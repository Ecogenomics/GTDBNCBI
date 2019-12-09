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

__prog_name__ = 'metadata_add_checkm_ncbi_to_database.py'
__prog_desc__ = 'Add CheckM information for NCBI genomes to database.'

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
import argparse
import tempfile
from collections import defaultdict

from biolib.taxonomy import Taxonomy


class AddCheckM(object):
    """Add CheckM information for NCBI genomes to database."""

    def __init__(self):
        self.metadata = {'Completeness': ['checkm_completeness', 'FLOAT'],
                         'Contamination': ['checkm_contamination', 'FLOAT'],
                         'Strain heterogeneity': ['checkm_strain_heterogeneity', 'FLOAT'],
                         'Marker lineage': ['checkm_marker_lineage', 'TEXT'],
                         '# genomes': ['checkm_genome_count', 'INT'],
                         '# markers': ['checkm_marker_count', 'INT'],
                         '# marker sets': ['checkm_marker_set_count', 'INT']}

    def run(self, checkm_profile_file, checkm_qa_sh100_file, genome_list_file):
        """Add CheckM data to database."""

        # get genomes to process
        genome_list = set()
        if genome_list_file:
            for line in open(genome_list_file):
                if line[0] == '#':
                    continue
                    
                if len(line.split('\t')) >= len(line.split(',')):
                    genome_list.add(line.rstrip().split('\t')[0])
                else:
                    genome_list.add(line.rstrip().split(',')[0])

        # add CheckM profile fields
        for header, data in self.metadata.iteritems():
            db_header, data_type = data

            num_genomes = 0
            temp_file = tempfile.NamedTemporaryFile(delete=False)
            with open(checkm_profile_file) as f:
                headers = f.readline().rstrip().split('\t')
                col_index = headers.index(header)

                for line in f:
                    line_split = line.rstrip().split('\t')
                    orig_genome_id = line_split[0]
                    genome_id = orig_genome_id[0:orig_genome_id.find('_', 4)]

                    if genome_id.startswith('GCA_'):
                        genome_id = 'GB_' + genome_id
                    elif genome_id.startswith('GCF_'):
                        genome_id = 'RS_' + genome_id

                    if genome_id not in genome_list:
                        print('Skipping genome: {}'.format(genome_id))
                        continue

                    data = line_split[col_index]
                    temp_file.write('%s\t%s\n' % (genome_id, data))
                    num_genomes += 1
                    
            if num_genomes == 0:
                print('No genomes identified.')
                sys.exit(-1)

            temp_file.close()
            cmd = 'gtdb -r metadata import --table %s --field %s --type %s --metadatafile %s' % ('metadata_genes', db_header, data_type, temp_file.name)
            print cmd
            os.system(cmd)
            os.remove(temp_file.name)
            
        # add strain heterogeneity results at 100%
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        with open(checkm_qa_sh100_file) as f:
            headers = f.readline().rstrip().split('\t')
            sh_index = headers.index('Strain heterogeneity')

            for line in f:
                line_split = line.rstrip().split('\t')
                genome_id = line_split[0]
                genome_id = genome_id[0:genome_id.find('_', 4)]

                if genome_id.startswith('GCA_'):
                    genome_id = 'GB_' + genome_id
                elif genome_id.startswith('GCF_'):
                    genome_id = 'RS_' + genome_id

                if genome_id not in genome_list:
                    continue

                sh = line_split[sh_index]
                temp_file.write('%s\t%s\n' % (genome_id, sh))

        temp_file.close()
        
        db_header = 'checkm_strain_heterogeneity_100'
        data_type = 'FLOAT'
        cmd = 'gtdb -r metadata import --table %s --field %s --type %s --metadatafile %s' % ('metadata_genes', db_header, data_type, temp_file.name)
        print cmd
        os.system(cmd)
        os.remove(temp_file.name)

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('checkm_profile_file', help='CheckM profile file for all genomes of interest')
    parser.add_argument('checkm_qa_sh100_file', help='CheckM QA file for 100%% strain heterogeneity for all genomes of interest')
    parser.add_argument('--genome_list', help='only process genomes in this list', default=None)

    args = parser.parse_args()

    try:
        p = AddCheckM()
        p.run(args.checkm_profile_file, args.checkm_qa_sh100_file, args.genome_list)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
