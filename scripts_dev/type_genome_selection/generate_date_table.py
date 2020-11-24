#!/usr/bin/env python3

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

__prog_name__ = 'generate_date_table.py'
__prog_desc__ = 'Generate table with LPSN year or priority for species and subspecies names.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2018'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

import os
import sys
import csv
import argparse
import re
import datetime
import logging

from biolib.logger import logger_setup


class DateEditor(object):
    """Main class
      """

    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')

    def parse_lpsn_scraped_priorities(self, lpsn_scraped_species_info):
        """Parse year of priority from references scraped from LPSN."""

        priorities = {}
        with open(lpsn_scraped_species_info) as lsi:
            lsi.readline()
            for line in lsi:
                infos = line.rstrip('\n').split('\t')
                
                sp = infos[0]
                if sp == 's__':
                    # *** hack to skip bad case in file
                    # Pierre to fix
                    continue 

                species_authority = infos[2]
                reference_str = species_authority.split(', ')[0]
                references = reference_str.replace('(', '').replace(')', '')
                years = re.sub(r'emend\.[^\d]*\d{4}', '', references)
                years = re.sub(r'ex [^\d]*\d{4}', ' ', years)
                years = re.findall('[1-3][0-9]{3}', years, re.DOTALL)
                years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]
                
                if len(years) == 0:
                    # assume this name is validated through ICN and just take the first 
                    # date given as the year of priority
                    years = re.findall('[1-3][0-9]{3}', references, re.DOTALL)
                    years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]
                    
                priorities[sp.replace('s__', '')] = years[0]
                            
        # We make sure that species and subspecies type species have the same date
        # ie Photorhabdus luminescens and Photorhabdus luminescens subsp.
        # Luminescens
        for k, v in priorities.items():
            infos_name = k.split(' ')
            if len(infos_name) == 2 and '{0} {1} subsp. {1}'.format(infos_name[0], infos_name[1]) in priorities:
                priorities[k] = min(int(v), int(priorities.get(
                    '{0} {1} subsp. {1}'.format(infos_name[0], infos_name[1]))))
            elif len(infos_name) == 4 and infos_name[1] == infos_name[3] and '{} {}'.format(infos_name[0], infos_name[1]) in priorities:
                priorities[k] = min(int(v), int(priorities.get(
                    '{} {}'.format(infos_name[0], infos_name[1]))))
                    
        return priorities
        
    def parse_lpsn_gss_priorities(self, lpsn_gss_file):
        """Get priority of species and usbspecies from LPSN GSS file."""

        priorities = {}
        illegitimate_names = set()
        with open(lpsn_gss_file, encoding='utf-8', errors='ignore') as f:
            csv_reader = csv.reader(f)

            for line_num, tokens in enumerate(csv_reader):
                if line_num == 0:
                    genus_idx = tokens.index('genus_name')
                    specific_idx = tokens.index('sp_epithet')
                    subsp_idx = tokens.index('subsp_epithet')
                    status_idx = tokens.index('status')
                    author_idx = tokens.index('authors')
                else:
                    generic = tokens[genus_idx].strip().replace('"', '')
                    specific = tokens[specific_idx].strip().replace('"', '')
                    subsp = tokens[subsp_idx].strip().replace('"', '')
                    
                    if subsp:
                        taxon = '{} {} subsp. {}'.format(generic, specific, subsp)
                    elif specific:
                        taxon = '{} {}'.format(generic, specific)
                    else:
                        # skip genus entries
                        continue

                    status = tokens[status_idx].strip().replace('"', '')
                    status_tokens = [t.strip() for t in status.split(';')]
                    status_tokens = [tt.strip() for t in status_tokens for tt in t.split(',') ]
                    
                    if 'illegitimate name' in status_tokens:
                        illegitimate_names.add(taxon)
                        if taxon in priorities:
                            continue

                    # get priority references, ignoring references if they are
                    # marked as being a revied name as indicated by a 'ex' or 'emend'
                    # (e.g. Holospora (ex Hafkine 1890) Gromov and Ossipov 1981)
                    ref_str = tokens[author_idx]
                    references = ref_str.replace('(', '').replace(')', '')
                    years = re.sub(r'emend\.[^\d]*\d{4}', '', references)
                    years = re.sub(r'ex [^\d]*\d{4}', ' ', years)
                    years = re.findall('[1-3][0-9]{3}', years, re.DOTALL)
                    years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]

                    if (taxon not in illegitimate_names
                        and taxon in priorities 
                        and years[0] != priorities[taxon]):
                            # conflict that can't be attributed to one of the entries being
                            # considered an illegitimate name
                            self.logger.error('Conflicting priority references for {}: {} {}'.format(
                                                taxon, years, priorities[taxon]))

                    priorities[taxon] = years[0]
        
        return priorities

    def run(self, lpsn_scraped_species_info, lpsn_gss_file, out_dir):
        """Parse priority year from LPSN data."""
        
        self.logger.info('Reading priority references scrapped from LPSN.')
        scraped_sp_priority = self.parse_lpsn_scraped_priorities(lpsn_scraped_species_info)
        self.logger.info(' - read priority for {:,} species.'.format(len(scraped_sp_priority)))
        
        self.logger.info('Reading priority references from LPSN GSS file.')
        gss_sp_priority = self.parse_lpsn_gss_priorities(lpsn_gss_file)
        self.logger.info(' - read priority for {:,} species.'.format(len(gss_sp_priority)))
        
        self.logger.info('Scrapped priority information for {:,} species not in GSS file.'.format(
                            len(set(scraped_sp_priority) - set(gss_sp_priority))))
        self.logger.info('Parsed priority information for {:,} species not on LPSN website.'.format(
                            len(set(gss_sp_priority) - set(scraped_sp_priority))))
                            
        self.logger.info('Writing out year of priority for species giving preference to GSS file.')
        output_file = open(os.path.join(out_dir, 'year_table.tsv'), 'w')
        same_year = 0
        diff_year = 0
        for sp in sorted(set(scraped_sp_priority).union(gss_sp_priority)):
            if sp in gss_sp_priority:
                output_file.write('{}\t{}\n'.format(sp, gss_sp_priority[sp]))
            else:
                output_file.write('{}\t{}\n'.format(sp, scraped_sp_priority[sp]))
                
            if sp in gss_sp_priority and sp in scraped_sp_priority:
                if gss_sp_priority[sp] == scraped_sp_priority[sp]:
                    same_year += 1
                else:
                    diff_year += 1
                    
        self.logger.info(' - same priority year in GSS file and website: {:,}'.format(same_year))
        self.logger.info(' - different priority year in GSS file and website: {:,}'.format(diff_year))
            
        output_file.close()


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--lpsn_scraped_species_info', 
                        help='LPSN species file created by LPSN website parsing.',
                        required=True)
    parser.add_argument('--lpsn_gss_file', 
                        help="table from lpsn.dsmz.de with nomenclature information (lpsn_gss_<date>.csv)",
                        required=True)
    parser.add_argument('--out_dir',
                        help='Output directory.',
                        required=True)

    args = parser.parse_args()
    
    logger_setup(args.out_dir, 
                __prog_name__.replace('.py', '.log'), 
                __prog_name__, 
                __version__, 
                False)
                
    try:
        dateeditor = DateEditor()
        dateeditor.run(args.lpsn_scraped_species_info,
                        args.lpsn_gss_file,
                        args.out_dir)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
        raise
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
