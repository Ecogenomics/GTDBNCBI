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

__prog_name__ = 'create_ssu_lsu_files.py'
__prog_desc__ = 'Create FASTA files with all 16S and 23S rRNA sequences from GTDB genomes.'

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
import gzip
import logging
import ntpath
import argparse
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy


class Create(object):
    """Create FASTA files with all 16S and 23S rRNA sequences from GTDB genomes."""

    def __init__(self):
        """Initialize."""

        pass

    def run(self, gtdb_bac_taxonomy_file, gtdb_ar_taxonomy_file, gtdb_path_file, gtdb_metadata_file, output_dir):
        """Create FASTA files with all 16S and 23S rRNA sequences from GTDB genomes."""
        
        # get User ID to UBA translation
        print('Reading GTDB metadata to translate User IDs to UBA IDs.')
        user_id_to_uba = {}
        with open(gtdb_metadata_file) as f:
            f.readline()
            
            for line in f:
                line_split = line.strip().split('\t')
                gid = line_split[0]
                org_name = line_split[1]
                if '(UBA' in org_name:
                    uba_id = org_name.split('(')[-1].replace(')', '')
                    user_id_to_uba[gid] = uba_id

        # read GTDB taxonomy
        print('Reading GTDB taxonomy.')
        gtdb_bac_taxonomy = Taxonomy().read(gtdb_bac_taxonomy_file)
        gtdb_ar_taxonomy = Taxonomy().read(gtdb_ar_taxonomy_file)
        gtdb_taxonomy = gtdb_bac_taxonomy.copy()
        gtdb_taxonomy.update(gtdb_ar_taxonomy)
                
        print('Identified %d bacterial genomes to process.' % len(gtdb_bac_taxonomy))
        print('Identified %d archaeal genomes to process.' % len(gtdb_ar_taxonomy))
        print('Identified %d genomes to process.' % len(gtdb_taxonomy))
        
        # read genome paths
        print('Reading path to genomes.')
        genome_paths = {}
        for line in open(gtdb_path_file):
            gid, gid_path = line.strip().split('\t')
            if gid in user_id_to_uba:
                gid = user_id_to_uba[gid]

            genome_paths[gid] = gid_path
                
        # sanity check data
        missing_paths = set(gtdb_taxonomy.keys()) - set(genome_paths.keys())
        if len(missing_paths) > 0:
            print('[WARNING] There are %d genomes in the taxonomy file without a specified genome path.' % len(missing_paths))

        # create FASTA file with 16S and 23S rRNA sequence files
        print('Parsing 16S and 23S rRNA sequence files.')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        fout_16S = open(os.path.join(output_dir, 'ssu.fna'), 'w')
        fout_23S = open(os.path.join(output_dir, 'lsu.fna'), 'w')
        missing_ssu = 0
        missing_lsu = 0
        for i, gid in enumerate(gtdb_taxonomy):
            if i % 1000 == 0:
                print('Processed %d genomes.' % i)
                
            if gid not in genome_paths:
                print('[WARNING] Genome %s does not have a specified genome path.' % gid)
                continue

            genome_path = genome_paths[gid]
            
            ssu_file = os.path.join(genome_path, 'rna_silva', 'ssu.fna')
            if not os.path.exists(ssu_file):
                missing_ssu += 1
                continue
                
            ssu_info_file = os.path.join(genome_path, 'rna_silva', 'ssu.hmm_summary.tsv')
            ssu_info = {}
            with open(ssu_info_file) as f:
                header = f.readline().strip().split('\t')
                contig_len_index = header.index('Sequence length')
                
                for line in f:
                    line_split = line.strip().split('\t')
                    
                    gene_id = line_split[0]
                    contig_length = int(line_split[contig_len_index])
                    ssu_info[gene_id] = contig_length
            
            for ssu_index, (seq_id, seq) in enumerate(seq_io.read_seq(ssu_file)):
                fout_16S.write('>%s~%s [ssu=%d bp] [contig=%d bp]\n' % (gid, seq_id, len(seq), ssu_info[seq_id]))
                fout_16S.write('%s\n' % seq)
                
            lsu_file = os.path.join(genome_path, 'rna_silva', 'lsu_23S.fna')
            if not os.path.exists(lsu_file):
                missing_lsu += 1
                continue
                
            lsu_info_file = os.path.join(genome_path, 'rna_silva', 'lsu_23S.hmm_summary.tsv')
            lsu_info = {}
            with open(lsu_info_file) as f:
                header = f.readline().strip().split('\t')
                contig_len_index = header.index('Sequence length')
                
                for line in f:
                    line_split = line.strip().split('\t')
                    
                    gene_id = line_split[0]
                    contig_length = int(line_split[contig_len_index])
                    lsu_info[gene_id] = contig_length
            
            for lsu_index, (seq_id, seq) in enumerate(seq_io.read_seq(lsu_file)):
                fout_23S.write('>%s~%s [ssu=%d bp] [contig=%d bp]\n' % (gid, seq_id, len(seq), lsu_info[seq_id]))
                fout_23S.write('%s\n' % seq)

        fout_16S.close()
        fout_23S.close()
        
        print('There were %d of %d (%.2f%%) genomes without an identifier 16S rRNA gene.' % (missing_ssu, 
                                                                                            len(gtdb_taxonomy),
                                                                                            missing_ssu*100.0/len(gtdb_taxonomy)))
                                                                                            
        print('There were %d of %d (%.2f%%) genomes without an identifier 23S rRNA gene.' % (missing_lsu, 
                                                                                            len(gtdb_taxonomy),
                                                                                            missing_lsu*100.0/len(gtdb_taxonomy)))
            
if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_bac_taxonomy_file', help='file with bacterial GTDB taxonomy for genomes to consider (download from GTDB website)')
    parser.add_argument('gtdb_ar_taxonomy_file', help='file with archaeal GTDB taxonomy for genomes to consider (download from GTDB website)')
    parser.add_argument('gtdb_path_file', help='file indicating path to GTDB genome data (dump from GTDB database)')
    parser.add_argument('gtdb_metadata_file', help='file indicating metadata for GTDB genome used to translate UBA identifiers (dump from GTDB database)')
    parser.add_argument('output_dir', help="output directory")

    args = parser.parse_args()

    try:
        p = Create()
        p.run(args.gtdb_bac_taxonomy_file, args.gtdb_ar_taxonomy_file, args.gtdb_path_file, args.gtdb_metadata_file, args.output_dir)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
