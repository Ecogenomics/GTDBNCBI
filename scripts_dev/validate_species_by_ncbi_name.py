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

__prog_name__ = 'validate_species_by_name.py'
__prog_desc__ = 'Compare GTDB species to NCBI species.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2017'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import re
import operator
import argparse
import tempfile
import shutil
from collections import defaultdict

from numpy import mean, std

from biolib.parallel import Parallel


class ValidateSpeciesNames(object):
    """Find discrepancies between GTDB and NCBI species names."""

    def __init__(self):
        pass

    def run(self, gtdb_taxonomy_file, ncbi_taxonomy_file):
        """Find discrepancies between GTDB and NCBI species names.."""
        
        # get NCBI species designations
        ncbi_species = {}
        for line in open(ncbi_taxonomy_file):
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            taxonomy = [taxon.strip() for taxon in line_split[1].split(';')]
            sp = taxonomy[6]
            ncbi_species[gid] = sp.replace('Candidatus ', '').replace('[', '').replace(']', '')

        # get GTDB species designations
        gtdb_species = {}
        unique_gtdb_species = set()
        for line in open(gtdb_taxonomy_file):
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            taxonomy = [taxon.strip() for taxon in line_split[1].split(';')]
            sp = taxonomy[6]

            # get canonical species names by removing polyphyly suffices
            canonical_sp = 's__' + ' '.join([t.strip() for t in re.split('_[A-Z]+', sp[3:])])
            gtdb_species[gid] = (sp, canonical_sp.strip())
            
            unique_gtdb_species.add(sp)
                
        print('Identified %d GTDB genomes.' % len(gtdb_species))
        print('There are %d GTDB species.' % len(unique_gtdb_species))
        
        discrepancies = {}
        bad_ncbi_sp = set()
        for gid in gtdb_species:
            gtdb_sp, gtdb_canonical_sp = gtdb_species[gid]
            ncbi_sp = ncbi_species[gid]
            
            if gtdb_canonical_sp != ncbi_sp and ncbi_sp != 's__':
                discrepancies[gid] = (gtdb_sp, gtdb_canonical_sp, ncbi_sp)
                bad_ncbi_sp.add(ncbi_sp)

        print('No. discrepancies', len(discrepancies))
        print('Bad NCBI species', len(bad_ncbi_sp))
                
        fout = open('gtdb_r80_species_discrepancies.tsv', 'w')
        fout.write('Genome ID\tGTDB species\tGTDB canonical\tNCBI species\n')
        
        sorted_d = sorted(discrepancies.items(), key=lambda x: x[1][1])
        for gid, info in sorted_d:
            fout.write('%s\t%s\n' % (gid, '\t'.join(info)))
        fout.close()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_taxonomy_file', help='file specifying GTDB taxonomy of all genomes')
    parser.add_argument('ncbi_taxonomy_file', help='file specifying NCBI taxonomy of all genomes')
  
    args = parser.parse_args()

    try:
        p = ValidateSpeciesNames()
        p.run(args.gtdb_taxonomy_file, 
                args.ncbi_taxonomy_file)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
