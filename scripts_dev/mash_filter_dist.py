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

__prog_name__ = 'mash_filter_dist.py'
__prog_desc__ = 'Filter pairwise Mash distances to reduce file size.'

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
import shutil
import tempfile
import argparse
from collections import defaultdict

from numpy import (mean as np_mean)

class MashFilter(object):
    """Filter pairwise Mash distances to reduce file size."""

    def __init__(self):
        """Initialization."""        
        pass

    def run(self, input_mash_dist_file, output_mash_dist_file, max_dist):
        """Filter pairwise Mash distances to reduce file size."""
	
	dists = defaultdict(list)
	for line in open(input_mash_dist_file):
		line_split = line.strip().split('\t')
		d = float(line_split[2])
		if d <= max_dist:
			gA = line_split[0]
			gB = line_split[1]
			if gA > gB:
				gA, gB = gB, gA
			dists[gA + '\t' + gB].append(d)
	
	fout = open(output_mash_dist_file, 'w')
	for k, d in dists.iteritems():
		fout.write('%s\t%.3f\n' % (k, np_mean(d)))
	fout.close()
	
if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_mash_dist_file', help='file with Mash distances to filter')
    parser.add_argument('output_mash_dist_file', help='file with Mash distances to filter')
    parser.add_argument('-d', '--max_dist', type=float, help='maximum distance to report', default=0.06)

    args = parser.parse_args()

    try:
        p = MashFilter()
        p.run(args.input_mash_dist_file, 
                args.output_mash_dist_file, 
                args.max_dist)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
