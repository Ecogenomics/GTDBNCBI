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
__prog_desc__ = 'Generate a tab delimited file listing all reference dates from LPSN,DSMZ,Straininfo for each species name.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2018'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

import sys
import argparse
import pickle
import os
import re
import datetime

# ('Terrimonas ferruginea', 'J Bacteriol 28, 415-431, 1934/Gwf-wasser-abwasser 117, 80-86, 1976', 1934)


class DateEditor(object):
    """Main class
      """

    def __init__(self):
        """Initialization."""
        
        pass

    def load_date_dict(self, species_info):
        """Parse year with priority from list of references."""
        
        if species_info is None:
            return {}
            
        datedict = {}
        with open(species_info) as lsi:
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
                    
                datedict[sp.replace('s__', '')] = years[0]
                            
        # We make sure that species and subspecies type species have the same date
        # ie Photorhabdus luminescens and Photorhabdus luminescens subsp.
        # Luminescens
        for k, v in datedict.items():
            infos_name = k.split(' ')
            if len(infos_name) == 2 and '{0} {1} subsp. {1}'.format(infos_name[0], infos_name[1]) in datedict:
                datedict[k] = min(int(v), int(datedict.get(
                    '{0} {1} subsp. {1}'.format(infos_name[0], infos_name[1]))))
            elif len(infos_name) == 4 and infos_name[1] == infos_name[3] and '{} {}'.format(infos_name[0], infos_name[1]) in datedict:
                datedict[k] = min(int(v), int(datedict.get(
                    '{} {}'.format(infos_name[0], infos_name[1]))))
                    
        return datedict

    def run(self, lpsn_species_info, outfile):
        """Parse priority year from LPSN data."""
        
        self.lpsn_date_dict = self.load_date_dict(lpsn_species_info)
        
        output_file = open(outfile, 'w')
        
        for spe in self.lpsn_date_dict:
            list_date = [''] * 4
            list_date[0] = spe
            if spe in self.lpsn_date_dict:
                list_date[1] = str(self.lpsn_date_dict.get(spe))

            output_file.write('{}\n'.format('\t'.join(list_date)))
            
        output_file.close()


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--lpsn_species_info', 
                        help='LPSN species file created by LPSN website parsing.',
                        required=True)
    parser.add_argument('--output',
                        help='Output file.',
                        required=True)

    args = parser.parse_args()

    try:
        dateeditor = DateEditor()
        dateeditor.run(args.lpsn_species_info,
                        args.output)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
        raise
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
