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

__prog_name__ = 'track_representatives.py'
__prog_desc__ = 'Identify changes to representative genomes between NCBI releases.'

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

csv.field_size_limit(sys.maxsize)

class TrackReps(object):
  """Identify changes to representative genomes between NCBI releases."""

  def __init__(self):
    pass

  def run(self, gtdb_metadata_prev, gtdb_metadata_cur):
    """Identify changes to representative genomes between NCBI releases."""
    

    # identify representative genomes in previous NCBI release
    print('Identifying representatives in previous NCBI release.')
    prev_reps = set()
    prev_genome_ids = set()
    header = True
    for row in csv.reader(open(gtdb_metadata_prev, 'rb')):
        if header:
            header = False
            
            gtdb_rep_index = row.index('gtdb_representative')
            gtdb_taxonomy_index = row.index('gtdb_taxonomy')
        else:
            genome_id = row[0]
            gtdb_taxonomy = row[gtdb_taxonomy_index]
            
            prev_genome_ids.add(genome_id)
            
            gtdb_rep = (row[gtdb_rep_index] == 't')
            if gtdb_rep:
                prev_reps.add(genome_id)
                
    print('%d representatives in previous NCBI release.' % len(prev_reps))
    
    # identify previous representatives in new NCBI release
    print('')
    print('Identifying representatives present in current NCBI release.')
    header = True
    cur_reps = set()
    cur_assigned_rep = {}
    current_genome_ids = set()
    uba_genome_ids = set()
    cur_genome_qual = {}
    fout = open('reps.tsv', 'w')
    for row in csv.reader(open(gtdb_metadata_cur, 'rb')):
        if header:
            header = False
            
            gtdb_rep_index = row.index('gtdb_representative')
            gtdb_genome_rep_index = row.index('gtdb_genome_representative')
            comp_index = row.index('checkm_completeness')
            cont_index = row.index('checkm_contamination')
        else:
            genome_id = row[0]
            organism_name = row[1]
            current_genome_ids.add(genome_id)
            
            rep_id = row[gtdb_genome_rep_index]
            cur_assigned_rep[genome_id] = rep_id
            
            comp = float(row[comp_index])
            cont = float(row[cont_index])
            cur_genome_qual[genome_id] = [comp, cont, comp-5*cont]
            
            if '(UBA' in organism_name:
                uba_genome_ids.add(genome_id)
            
            gtdb_rep = (row[gtdb_rep_index] == 't')
            if gtdb_rep:
                cur_reps.add(genome_id)
                
            if genome_id in prev_reps:
                fout.write('%s\t%s\n' % (genome_id, genome_id))
                
    fout.close()
                
    print('%d representatives in current NCBI release.' % len(cur_reps))
    
    print('%d UBA genomes.' % len(uba_genome_ids))
    
    print('')
    print('%d representatives in common between NCBI releases.' % len(cur_reps.intersection(prev_reps)))
    missing_reps = prev_reps - current_genome_ids
    print('%d previous representatives no longer in NCBI.' % len(missing_reps))
    deprecated_reps = prev_reps.intersection(current_genome_ids) - cur_reps
    print('%d previous representatives no longer a representative.' % len(deprecated_reps))
    
    user_genome_rep = set()
    for genome_id in deprecated_reps:
        if genome_id.startswith('U_') and genome_id not in uba_genome_ids:
            user_genome_rep.add(genome_id)
            
    print('  %d of these are deprecated User genomes.' % len(user_genome_rep))
    
    for genome_id in (deprecated_reps - user_genome_rep):
        comp, cont, q = cur_genome_qual[genome_id]
        rep_id = cur_assigned_rep[genome_id]
        rep_comp, rep_cont, rep_q = cur_genome_qual[rep_id]
        print(genome_id, comp, cont, q, rep_id, rep_comp, rep_cont, rep_q)
    
    print('')
    elevated_to_reps = cur_reps.intersection(prev_genome_ids) - prev_reps
    print('%d genomes elevated to a representatives.' % len(elevated_to_reps))
    new_reps = cur_reps - prev_genome_ids
    print('%d representatives absent from previous release.' % len(new_reps))
    
    #print '%d representatives absent from current NCBI release.' % len(remaining_prev_reps)
    
    # try to identify what happened to absent representatives
    if False:
        print('')
        print('Identifying status of absent representatives.')
        
        remaining_prev_reps = prev_reps - cur_reps
        
        moved_to_refseq = set()
        moved_to_genbank = set()
        new_genome_version = set()
        for genome_id in current_genome_ids:
            if genome_id.startswith('U_'):
                continue
                
            # check for database or version change
            cur_version = int(genome_id.split('.')[-1])
            for new_version in range(1, cur_version+5):
                new_version_id = genome_id.replace('.%d' % cur_version, '.%d' % new_version)
                if new_version_id in remaining_prev_reps:
                    new_genome_version.add(new_version_id)
                    continue
            
                gb_genome_id = new_version_id.replace('RS_GCF', 'GB_GCA')
                if gb_genome_id in remaining_prev_reps:
                    moved_to_refseq.add(gb_genome_id)
                    continue
                
                rs_genome_id = new_version_id.replace('GB_GCA', 'RS_GCF')
                if rs_genome_id in remaining_prev_reps:
                    moved_to_genbank.add(rs_genome_id)
                    continue

        print('%d representatives moved from GenBank to RefSeq.' % len(moved_to_genbank))
        print('%d representatives moved from RefSeq to GenBank.' % len(moved_to_refseq))
        print('%d representatives have an updated version number.' % len(new_genome_version))
        
        remaining_prev_reps = remaining_prev_reps - moved_to_genbank - moved_to_refseq - new_genome_version
        print('Remaining %d representatives:' % len(remaining_prev_reps))
        #print remaining_prev_reps

if __name__ == '__main__':
  print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
  print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('gtdb_metadata_prev', help='GTDB metadata for previous NCBI release.')
  parser.add_argument('gtdb_metadata_cur', help='GTDB metadata for current NCBI release.')

  args = parser.parse_args()

  try:
    p = TrackReps()
    p.run(args.gtdb_metadata_prev, args.gtdb_metadata_cur)
  except SystemExit:
    print("\nControlled exit resulting from an unrecoverable error or warning.")
  except:
    print("\nUnexpected error:", sys.exc_info()[0])
    raise
