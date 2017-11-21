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

__prog_name__ = 'validate_type_species.py'
__prog_desc__ = 'Check for LPSN and DSMZ type species missing in GTDB taxonomy.'

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
import argparse
import tempfile
import shutil
from collections import defaultdict, Counter


class ValidateTypeSpecies(object):
    """Check for LPSN and DSMZ type species missing in GTDB taxonomy."""

    def __init__(self):
        pass

    def run(self, gtdb_metadata_file, lpsn_species, output_dir):
        """Check for LPSN and DSMZ type species missing in GTDB taxonomy."""
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        # parse LPSN type species information
        lpsn_type_species = set()
        with open(lpsn_species) as f:
            f.readline()
            
            for line in f:
                line_split = line.strip().split('\t')
                
                lpsn_sp = line_split[0]
                lpsn_type_sp = line_split[1]
                lpsn_authority = line_split[2]
                
                if lpsn_type_sp:
                    if lpsn_sp in lpsn_type_species:
                        print lpsn_sp
                    lpsn_type_species.add(lpsn_sp)
                    
        print 'Identified %d type species at LPSN.' % len(lpsn_type_species)

        # get NCBI species designations
        gtdb_species = set()
        lpsn_type_sp_at_ncbi = set()
        lpsn_type_sp_at_gtdb = set()
        type_species_under_gtdb = defaultdict(set)
        gtdb_type_species = defaultdict(set)
        with open(gtdb_metadata_file) as f:
            header = f.readline().strip().split('\t')
            
            gtdb_taxonomy_index = header.index('gtdb_taxonomy')
            ncbi_taxonomy_index = header.index('ncbi_taxonomy')

            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[0]
                
                gtdb_taxonomy = line_split[gtdb_taxonomy_index]
                gtdb_sp = None
                if gtdb_taxonomy:
                    gtdb_taxa = [t.strip() for t in gtdb_taxonomy.split(';')]
                    gtdb_sp = gtdb_taxa[6]
                    #canonical_gtdb_sp = 's__' + ' '.join([t.strip() for t in re.split('_[A-Z]+', gtdb_sp[3:])]).strip()
                    gtdb_species.add(gtdb_sp)

                ncbi_taxonomy = line_split[ncbi_taxonomy_index]
                ncbi_sp = None
                if ncbi_taxonomy and ncbi_taxonomy != 'none':
                    ncbi_taxa = [t.strip() for t in ncbi_taxonomy.split(';')]
                    ncbi_sp = ncbi_taxa[6]
                    
                if gid == 'GB_GCA_900101375.1':
                    print 'GB_GCA_900101375.1', ncbi_sp
                    
                # If NCBI has a genome marked as being from the type species,
                # then we would expect there to be a GTDB genome assigned to this type species.
                # We would also expect this to be a genome assigned to the type species at NCBI.
                if ncbi_sp in lpsn_type_species:
                    lpsn_type_sp_at_ncbi.add(ncbi_sp)
                    type_species_under_gtdb[ncbi_sp].add(str(gtdb_sp))
                        
                if gtdb_sp in lpsn_type_species:
                    lpsn_type_sp_at_gtdb.add(gtdb_sp)
                    gtdb_type_species[gtdb_sp].add(str(ncbi_sp))
                    
        print 'Identified %d LPSN type species at NCBI.' % len(lpsn_type_sp_at_ncbi)
        print 'Identified %d LPSN type species in GTDB.' % len(lpsn_type_sp_at_gtdb)
        print 'LPSN type species missing in GTDB: %s' % len(lpsn_type_sp_at_ncbi - lpsn_type_sp_at_gtdb)
        
        fout = open(os.path.join(output_dir, 'suspicious_lpsn_type_species.tsv'), 'w')
        fout.write('Type species\tMissing in GTDB\tAssignments under GTDB taxonomy\tAssignments under NCBI taxonomy\n')
        for lpsn_sp in lpsn_type_sp_at_ncbi - lpsn_type_sp_at_gtdb:
            fout.write('%s\t%s\t%s\t%s\n' % (lpsn_sp, 
                                                'True', 
                                                ','.join(type_species_under_gtdb.get(lpsn_sp, [])),
                                                'N/A'))

        for gtdb_sp, ncbi_sps in gtdb_type_species.iteritems():
            if gtdb_sp not in ncbi_sps:
                fout.write('%s\t%s\t%s\t%s\n' % (gtdb_sp, 'False', 'N/A', ','.join(ncbi_sps)))
        fout.close()
        
        
if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_metadata_file', help='file specifying GTDB metadata (TSV)')
    parser.add_argument('lpsn_species', help='file with LPSN type species information (TSV)')
    parser.add_argument('output_dir', help='output directory')
  
    args = parser.parse_args()

    try:
        p = ValidateTypeSpecies()
        p.run(args.gtdb_metadata_file,
                args.lpsn_species,
                args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
