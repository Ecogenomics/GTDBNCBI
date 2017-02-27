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

__prog_name__ = 'propagate_gtdb_taxonomy.py'
__prog_desc__ = 'Propagate GTDB taxonomy between NCBI releases.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2016'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import csv
import argparse
from collections import defaultdict

class Propagate(object):
  """Propagate GTDB taxonomy between NCBI releases."""

  def __init__(self):
    pass

  def run(self, gtdb_metadata_prev, gtdb_metadata_cur, taxonomy_file, new_rep_file):
    """Propagate GTDB taxonomy between NCBI releases."""
    
    # get GTDB taxonomy for genome in previous release
    print 'Reading GTDB taxonomy of genome in previous release:'
    prev_gtdb_taxonomy = {}
    prev_gtdb_genomes = set()
    prev_is_rep = set()
    header = True
    for row in csv.reader(open(gtdb_metadata_prev, 'rb')):
        if header:
            header = False
            
            gtdb_taxonomy_index = row.index('gtdb_taxonomy')
            gtdb_rep_index = row.index('gtdb_representative')
        else:
            genome_id = row[0]
            prev_gtdb_genomes.add(genome_id)
            
            gtdb_taxonomy = row[gtdb_taxonomy_index]
            if gtdb_taxonomy:
                prev_gtdb_taxonomy[genome_id] = gtdb_taxonomy
                
            is_rep = (row[gtdb_rep_index] == 't')
            if is_rep:
                prev_is_rep.add(genome_id)
                
    print '  %d of %d (%.1f%%) genomes in previous NCBI release had a GTDB taxonomy string' % (len(prev_gtdb_taxonomy), 
                                                                                        len(prev_gtdb_genomes),
                                                                                        len(prev_gtdb_taxonomy)*100.0/len(prev_gtdb_genomes))
                                                                                        
    print '  %d genomes were identified as representatives' % len(prev_is_rep)
    
    # identify previous representatives in new NCBI release
    print ''
    print 'Identifying unchanged genomes in current NCBI release:'
    header = True
    fout = open(taxonomy_file, 'w')
    retained_genomes = set()
    current_genome_ids = []
    prev_rep_count = 0
    for row in csv.reader(open(gtdb_metadata_cur, 'rb')):
        if header:
            header = False
            
            gtdb_rep_index = row.index('gtdb_representative')
        else:
            genome_id = row[0]
            current_genome_ids.append(genome_id)
            
            if genome_id in prev_gtdb_genomes:
                retained_genomes.add(genome_id)  
                if genome_id in prev_gtdb_taxonomy:
                    fout.write('%s\t%s\n' % (genome_id, prev_gtdb_taxonomy[genome_id]))
                    
                if genome_id in prev_is_rep:
                    prev_rep_count += 1

    remaining_prev_genomes = prev_gtdb_genomes - retained_genomes
    print '  %d (%.1f%%) genomes unchanged in current NCBI release' % (len(retained_genomes),
                                                                        len(retained_genomes)*100.0/len(prev_gtdb_genomes))
    print '  %d (%.1f%%) genomes absent or modified in current NCBI release' % (len(remaining_prev_genomes),
                                                                        len(remaining_prev_genomes)*100.0/len(prev_gtdb_genomes))
    print '  %d representatives unchanged in current GTDB release' % prev_rep_count
    
    # try to identify what happened to absent representatives
    print ''
    print 'Identifying genomes that have changed databases or version:'
    
    fout_new_reps = open(new_rep_file, 'w')
    moved_to_refseq = set()
    moved_to_genbank = set()
    new_genome_version = set()
    for genome_id in current_genome_ids:
        if genome_id.startswith('U_'):
            continue
            
        # check for database or version change
        cur_version = int(genome_id.split('.')[-1])
        for new_version in xrange(1, cur_version+5):
            new_version_id = genome_id.replace('.%d' % cur_version, '.%d' % new_version)
            if new_version_id in remaining_prev_genomes:
                new_genome_version.add(new_version_id)
                if new_version_id in prev_gtdb_taxonomy:
                    fout.write('%s\t%s\n' % (genome_id, prev_gtdb_taxonomy[new_version_id]))

                if new_version_id in prev_is_rep:
                    fout_new_reps.write('%s\t%s\n' % (genome_id, str(True)))
                
                continue
        
            gb_genome_id = new_version_id.replace('RS_GCF', 'GB_GCA')
            if gb_genome_id in remaining_prev_genomes:
                moved_to_refseq.add(gb_genome_id)
                if gb_genome_id in prev_gtdb_taxonomy:
                    fout.write('%s\t%s\n' % (genome_id, prev_gtdb_taxonomy[gb_genome_id]))
                
                if gb_genome_id in prev_is_rep:
                    fout_new_reps.write('%s\t%s\n' % (genome_id, str(True)))
                    
                continue
            
            rs_genome_id = new_version_id.replace('GB_GCA', 'RS_GCF')
            if rs_genome_id in remaining_prev_genomes:
                moved_to_genbank.add(rs_genome_id)
                if rs_genome_id in prev_gtdb_taxonomy:
                    fout.write('%s\t%s\n' % (genome_id, prev_gtdb_taxonomy[rs_genome_id]))
                
                if rs_genome_id in prev_is_rep:
                    fout_new_reps.write('%s\t%s\n' % (genome_id, str(True)))
                    
                continue

    print '  %d (%.1f%%) genomes moved from GenBank to RefSeq' % (len(moved_to_genbank), len(moved_to_genbank)*100.0/len(prev_gtdb_genomes))
    print '  %d (%.1f%%) genomes moved from RefSeq to GenBank' % (len(moved_to_refseq), len(moved_to_refseq)*100.0/len(prev_gtdb_genomes))
    print '  %d (%.1f%%) genomes have a new version number' % (len(new_genome_version), len(new_genome_version)*100.0/len(prev_gtdb_genomes))
    
    remaining_prev_genomes = remaining_prev_genomes - moved_to_genbank - moved_to_refseq - new_genome_version
    print ''
    print 'There are %d genomes not present in the current release.' % len(remaining_prev_genomes)
    print remaining_prev_genomes
    
    fout.close()
    fout_new_reps.close()
    
if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('gtdb_metadata_prev', help='GTDB metadata for previous NCBI release.')
  parser.add_argument('gtdb_metadata_cur', help='GTDB metadata for current NCBI release.')
  parser.add_argument('taxonomy_file', help='Propagated GTDB taxonomy file.')
  parser.add_argument('new_rep_file', help='GTDB representatives with modified accession IDs.')

  args = parser.parse_args()

  try:
    p = Propagate()
    p.run(args.gtdb_metadata_prev, args.gtdb_metadata_cur, args.taxonomy_file, args.new_rep_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
