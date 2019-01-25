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

    def __init__(self, lpsn_species_info, dsmz_species_info, straininfo_species_info):
        """Initialization."""
        self.year = datetime.datetime.now().year
        self.lpsn_date_dict = self.load_date_dict(lpsn_species_info)
        self.dsmz_date_dict = self.load_date_dict(dsmz_species_info)
        self.straininfo_date_dict = self.load_date_dict(
            straininfo_species_info)

    def load_date_dict(self, lpsn_species_info):
        datedict = {}
        with open(lpsn_species_info) as lsi:
            for line in lsi:
                if 'straininfo_strains_number' in line:
                    continue
                infos = line.rstrip('\n').split('\t')
                if infos[2] != '':
                    date = re.sub(r'\([^)]*\)', '', infos[2])
                    date = re.sub(r'emend\.[^\d]*\d{4}', '', date)
                    date = re.sub(r' ex [^\d]*\d{4}', ' ', date)
                    matches = re.findall('[1-3][0-9]{3}', date, re.DOTALL)
                    matches = [int(a) for a in matches if int(a) <= self.year]
                    if matches:
                        datedict[infos[0].replace('s__', '')] = sorted(
                            matches)[0]
        # We make sure that species and subspecies type species have the same date
        # ie Photorhabdus luminescens and Photorhabdus luminescens subsp.
        # Luminescens
        for k, v in datedict.iteritems():
            infos_name = k.split(' ')
            if len(infos_name) == 2 and '{0} {1} subsp. {1}'.format(infos_name[0], infos_name[1]) in datedict:
                datedict[k] = min(int(v), int(datedict.get(
                    '{0} {1} subsp. {1}'.format(infos_name[0], infos_name[1]))))
            elif len(infos_name) == 4 and infos_name[1] == infos_name[3] and '{} {}'.format(infos_name[0], infos_name[1]) in datedict:
                datedict[k] = min(int(v), int(datedict.get(
                    '{} {}'.format(infos_name[0], infos_name[1]))))
        return datedict

    def run(self, outfile):
        output_file = open(outfile, 'w')
        all_species = self.lpsn_date_dict.keys()
        all_species.extend(self.dsmz_date_dict.keys())
        all_species.extend(self.straininfo_date_dict.keys())
        for spe in set(all_species):
            list_date = [''] * 4
            list_date[0] = spe
            if spe in self.lpsn_date_dict:
                list_date[1] = str(self.lpsn_date_dict.get(spe))
            if spe in self.dsmz_date_dict:
                list_date[2] = str(self.dsmz_date_dict.get(spe))
            if spe in self.straininfo_date_dict:
                list_date[3] = str(self.straininfo_date_dict.get(spe))
            output_file.write('{}\n'.format('\t'.join(list_date)))
        output_file.close()


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--lpsn_species_info', help='LPSN species file created by LPSN website parsing.')
    parser.add_argument(
        '--dsmz_species_info', help='DSMZ species file created by DSMZ website parsing.')
    parser.add_argument(
        '--straininfo_species_info', help='Straininfo strain file created by Straininfo website parsing.')
    parser.add_argument('--output',
                        help='Output file.')

    args = parser.parse_args()

    try:
        dateeditor = DateEditor(
            args.lpsn_species_info, args.dsmz_species_info,
            args.straininfo_species_info)
        dateeditor.run(args.output)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
