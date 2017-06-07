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

__prog_name__ = 'set_gtdb_taxonomy_to_reps.py'
__prog_desc__ = 'Force clustered genomes to have same GTDB taxonomy.'

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
import argparse
import logging
from collections import defaultdict, Counter

from biolib.taxonomy import Taxonomy

from gtdb import GenomeDatabase
from gtdb.Exceptions import (GenomeDatabaseError, 
                                DumpDBErrors, 
                                DumpDBWarnings, 
                                ErrorReport)
                                
import psycopg2
from psycopg2.extensions import AsIs


class Script():
    def __init__(self):
        """Initialization."""
        
        logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S",
                            level=logging.DEBUG)
        self.logger = logging
        
    def setup_db(self, gtdb_version, threads):
        """Setup database."""
        
        # initialise the backend
        self.db = GenomeDatabase.GenomeDatabase(threads, False)
        self.db.conn.MakePostgresConnection(gtdb_version)

        # login
        try:
            self.db.Login(None, False)
        except GenomeDatabaseError as e:
            self.db.conn.ClosePostgresConnection()
            ErrorReport(e.message + " The following error(s) were reported:\n")
            DumpDBErrors(self.db)
            sys.exit(-1)
            
    def run(self, threads):
        """Run script."""
        
        self.logger.info('Determining GTDB representatives and their taxonomic information.')
        
        # get representatives and clustered genomes
        cur = self.db.conn.cursor()
        
        q = ("SELECT id, genome, gtdb_taxonomy, ncbi_taxonomy, gtdb_genome_representative FROM metadata_view")
        cur.execute(q)
        
        reps = {}
        clustered = {}
        clusters = defaultdict(set)
        gtdb_taxa = {}
        ncbi_taxa = {}
        reps_without_taxonomy = set()
        for gid, genome, gtdb_taxonomy, ncbi_taxonomy, gtdb_genome_representative in cur:
            if not gtdb_genome_representative:
                continue
                
            if ncbi_taxonomy:
                ncbi_taxa[genome] = map(str.strip, ncbi_taxonomy.split(';'))
                
            if gtdb_taxonomy:
                gtdb_taxa[genome] = map(str.strip, gtdb_taxonomy.split(';'))
            else:
                gtdb_taxa[genome] = list(Taxonomy.rank_prefixes)
                if ncbi_taxonomy:
                    # fill in domain
                    gtdb_taxa[genome][0] = ncbi_taxa[genome][0]
                
            if genome == gtdb_genome_representative:
                if not gtdb_taxonomy or gtdb_taxa[genome][1] == 'p__':
                    reps[genome] = ';'.join(gtdb_taxa[genome])
                    reps_without_taxonomy.add(genome)
                else:
                    reps[genome] = gtdb_taxonomy
                    
            else:
                clustered[gid] = gtdb_genome_representative
                
                clusters[gtdb_genome_representative].add(genome)
                   
        self.logger.info('Identified %d GTDB representatives.' % len(reps))
        self.logger.info('  %d reps have no assigned GTDB taxonomy' % len(reps_without_taxonomy))
        self.logger.info('Identified %d clustered genomes.' % len(clustered))
        
        # assign representatives without a taxonomy string the
        # consensus taxonomy of their cluster
        assigned_taxonomy = 0
        for rep_id in reps_without_taxonomy:
            if not len(clusters[rep_id]):
                continue
            
            taxa = []
            for r in xrange(0, 7):
                t = []
                for gid in clusters[rep_id]:
                    t.append(gtdb_taxa[gid][r])

                taxon, count = Counter(t).most_common(1)[0]
                taxa.append(taxon)
                
            if taxa[0] == 'd__' and rep_id in ncbi_taxa:
                taxa[0] = ncbi_taxa[rep_id][0]
                print rep_id, taxa[0]
            
            reps[rep_id] = ';'.join(taxa)
            assigned_taxonomy += 1
            
        self.logger.info('Assigned %d representatives a consensus taxonomy.' % assigned_taxonomy)
        
        clustered_tax = []
        for gid, rep_id in clustered.iteritems():
            gtdb_taxonomy = reps[rep_id]
            gtdb_taxa = map(str.strip, gtdb_taxonomy.split(';'))
            d = [gtdb_taxonomy] + gtdb_taxa + [gid]
            clustered_tax.append(d)
                
        self.logger.info('Assigning representative taxonomy to %d clustered genomes.' % len(clustered_tax))
            
        q = ('UPDATE metadata_taxonomy SET gtdb_taxonomy = %s, '
                + 'gtdb_domain = %s, gtdb_phylum = %s, '
                + 'gtdb_class = %s, gtdb_order = %s, '
                + 'gtdb_family = %s, gtdb_genus = %s, '
                + 'gtdb_species = %s '
                + 'WHERE id = %s')
        cur.executemany(q, clustered_tax)
            
        self.db.conn.commit()
        cur.close()
 
if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_version', help='GTDB database version (i.e., gtdb_releaseX)')
    parser.add_argument('-t', '--threads', help='number of CPUs to use', type=int, default=1)

    args = parser.parse_args() 

    try:
        p = Script()
        
        # setup database
        p.setup_db(args.gtdb_version, args.threads)
    
        p.run(args.threads)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise