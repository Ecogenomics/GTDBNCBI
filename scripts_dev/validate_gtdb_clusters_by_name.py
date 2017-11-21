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

__prog_name__ = 'validate_gtdb_clusters_by_name.py'
__prog_desc__ = 'Tabulate NCBI species assigned to each GTDB cluster.'

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


class ValidateClusters(object):
    """Tabulate NCBI species assigned to each GTDB cluster."""

    def __init__(self):
        pass

    def run(self, gtdb_metadata_file, output_dir):
        """Tabulate NCBI species assigned to each GTDB cluster."""
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        clusters = defaultdict(set)
        gtdb_sp = {}
        ncbi_sp = {}
        gtdb_taxonomy = {}
        with open(gtdb_metadata_file) as f:
            header = f.readline().strip().split('\t')
            
            gtdb_taxonomy_index = header.index('gtdb_taxonomy')
            ncbi_taxonomy_index = header.index('ncbi_taxonomy')
            rep_id_index = header.index('gtdb_genome_representative')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[0]
                rep_id = line_split[rep_id_index]
                if not rep_id or rep_id == 'none':
                    continue
                
                ncbi_taxonomy = line_split[ncbi_taxonomy_index]
                if ncbi_taxonomy and ncbi_taxonomy != 'none':
                    ncbi_taxa = [t.strip() for t in ncbi_taxonomy.split(';')]
                    ncbi_sp[gid] = ncbi_taxa[6]
                    
                tmp_gtdb_taxonomy = line_split[gtdb_taxonomy_index]
                if tmp_gtdb_taxonomy and tmp_gtdb_taxonomy != 'none':
                    gtdb_taxa = [t.strip() for t in tmp_gtdb_taxonomy.split(';')]
                    gtdb_sp[gid] = gtdb_taxa[6]
                    gtdb_taxonomy[gid] = tmp_gtdb_taxonomy
                    
                clusters[rep_id].add(gid)
                
        fout = open(os.path.join(output_dir, 'gtdb_clusters_info.tsv'), 'w')
        fout.write('Representative ID\tNo. genomes\tGTDB taxonomy\tGTDB species\tNCBI species')
        fout.write('\tNo. NCBI species\tNo. NCBI species >1\tNo. incongruent with GTDB name\n')
        for rep_id in sorted(clusters, key=lambda k: len(clusters[k]), reverse=True):
            fout.write('%s\t%d\t%s' % (rep_id, len(clusters[rep_id]), gtdb_taxonomy[rep_id]))
            gtdb_species = defaultdict(int)
            ncbi_species = defaultdict(int)
            for gid in clusters[rep_id]:
                if gid in gtdb_sp:
                    gtdb_species[gtdb_sp[gid]] += 1
                    
                if gid in ncbi_sp:
                    ncbi_species[ncbi_sp[gid]] += 1
                    
            if len(gtdb_species) != 1:
                print '[Error] Cluster assigned to more than one GTDB species: %s' % rep_id
                sys.exit()
                
            gtdb_key = gtdb_species.keys()[0]
            canonical_gtdb_sp = 's__' + ' '.join([t.strip() for t in re.split('_[A-Z]+', gtdb_key[3:])]).strip()
            fout.write('\t%s: %d' % (gtdb_key, gtdb_species[gtdb_key]))
                
            ncbi_str = []
            num_ncbi_species = 0
            num_ncbi_species_multi_genomes = 0
            num_incongruent = 0
            for sp in sorted(ncbi_species, key=lambda k: ncbi_species[k], reverse=True):
                ncbi_str.append('%s: %d' % (sp, ncbi_species[sp]))
                if sp != 's__':
                    num_ncbi_species += 1
                    if ncbi_species[sp] > 1:
                        num_ncbi_species_multi_genomes += 1
                        if sp != canonical_gtdb_sp:
                            num_incongruent += ncbi_species[sp]
            
            fout.write('\t%s\t%d\t%d\t%d\n' % (', '.join(ncbi_str), 
                                            num_ncbi_species, 
                                            num_ncbi_species_multi_genomes, 
                                            num_incongruent))
        fout.close()
            

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_metadata_file', help='file specifying GTDB metadata (TSV file)')
    parser.add_argument('output_dir', help='output directory')
    
    args = parser.parse_args()

    try:
        p = ValidateClusters()
        p.run(args.gtdb_metadata_file,
                args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
