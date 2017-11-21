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

import os
import sys
import argparse

__prog_name__ = 'straininfo_scraper.py'
__prog_desc__ = 'Produce metadata files describing type genera, species, and strains according to Straininfo'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'


class PNUClient(object):

    def run(self, infile, outdir):
        dict_species = {}
        with open(infile) as strainfile:
            strainfile.readline()
            for line in strainfile:
                info = line.strip().split(';')
                if len(info[2].split()) != 2:
                    continue
                if info[2] in dict_species:
                    dict_species.get(info[2]).append(info[1])
                else:
                    dict_species[info[2]] = [info[1]]

            outfile = open(os.path.join(outdir, 'straininfo_strains.txt'), 'w')
            outfile.write("straininfo_strains\n")
            for k, v in dict_species.iteritems():
                outfile.write("{0} {1}\n".format(k, "=".join(set(v))))


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('straininfo', help='straininfo file')
    parser.add_argument('output_dir', help='output directory')

    args = parser.parse_args()

    try:
        p = PNUClient()
        p.run(args.straininfo, args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
