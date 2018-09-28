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

__prog_name__ = 'silva_gtdb_assignments.py'
__prog_desc__ = 'Create table assigning GTDB taxonomy to SILVA accessions based on SSU and LSU BLAST results.'

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
    """Create table assigning GTDB taxonomy to SILVA accessions based on SSU and LSU BLAST results."""

    def __init__(self):
        """Initialize."""
        
        self.per_iden_threshold = 98.7
        self.per_aln_len_threshold = 98.7

        self.min_ssu_len = {'d__Bacteria': 1200, 'd__Archaea': 900}
        self.min_lsu_len = {'d__Bacteria': 1900, 'd__Archaea': 1900}
        
    def _parse_blast_table(self, ssu_blast_table, gtdb_taxonomy, silva_ssu_taxonomy, min_lens, output_table):
        """Parse BLAST results."""
        
        fout = open(output_table, 'w')
        fout.write('GTDB ID\tSILVA ID\tGTDB taxonomy\tSILVA taxonomy\tLength (bp)\tPercent identity\tPercent alignment length\tE-value\n')
        for line in open(ssu_blast_table):
            line_split = line.strip().split('\t')
            
            qid = line_split[0]
            gid = qid.split('~')[0]
            gtdb_domain = gtdb_taxonomy[gid][0]

            qlen = int(line_split[1])
            if qlen < min_lens[gtdb_domain]:
                continue
                
            sid = line_split[2]
            slen = int(line_split[3])
            
            aln_len = int(line_split[4])
            mismatch = int(line_split[5])
            gaps = int(line_split[6])
            
            pident = float(line_split[7])
            evalue = float(line_split[9])
            
            adj_aln_len = aln_len - gaps
            per_aln_len = adj_aln_len * 100.0 / qlen
            per_iden = (adj_aln_len - mismatch) * 100.0 / qlen

            if (per_iden >= self.per_iden_threshold
                    and per_aln_len >= self.per_aln_len_threshold):
                        fout.write('%s\t%s\t%s\t%s\t%d\t%.2f\t%.2f\t%g\n' % (qid, 
                                                                                sid, 
                                                                                '; '.join(gtdb_taxonomy[gid]),
                                                                                silva_ssu_taxonomy[sid],
                                                                                qlen,
                                                                                per_iden, 
                                                                                per_aln_len, 
                                                                                evalue))

        fout.close()

    def run(self, 
                gtdb_bac_taxonomy_file, 
                gtdb_ar_taxonomy_file, 
                silva_ssu_ref,
                silva_lsu_ref,
                ssu_blast_table, 
                lsu_blast_table, 
                output_dir):
        """Create table assigning GTDB taxonomy to SILVA accessions based on SSU and LSU BLAST results."""
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # read GTDB taxonomy
        print('Reading GTDB taxonomy.')
        gtdb_bac_taxonomy = Taxonomy().read(gtdb_bac_taxonomy_file)
        gtdb_ar_taxonomy = Taxonomy().read(gtdb_ar_taxonomy_file)
        gtdb_taxonomy = gtdb_bac_taxonomy.copy()
        gtdb_taxonomy.update(gtdb_ar_taxonomy)
        
        print('Identified %d bacterial genomes to process.' % len(gtdb_bac_taxonomy))
        print('Identified %d archaeal genomes to process.' % len(gtdb_ar_taxonomy))
        print('Identified %d genomes to process.' % len(gtdb_taxonomy))
        
        # read SILVA taxonomy
        print('Reading SILVA 16S and 23S rRNA taxonomies.')
        silva_ssu_taxonomy = {}
        for seq_id, seq, taxonomy in seq_io.read_seq(silva_ssu_ref, keep_annotation=True):
            silva_ssu_taxonomy[seq_id] = taxonomy
            
        silva_lsu_taxonomy = {}
        for seq_id, seq, taxonomy in seq_io.read_seq(silva_lsu_ref, keep_annotation=True):
            silva_lsu_taxonomy[seq_id] = taxonomy

        # parse BLAST tables
        print('Parsing BLAST tables.')
        
        ssu_table = os.path.join(output_dir, 'ssu_silva.tsv')
        self._parse_blast_table(ssu_blast_table, gtdb_taxonomy, silva_ssu_taxonomy, self.min_ssu_len, ssu_table)
        
        lsu_table = os.path.join(output_dir, 'lsu_silva.tsv')
        self._parse_blast_table(lsu_blast_table, gtdb_taxonomy, silva_lsu_taxonomy, self.min_lsu_len, lsu_table)

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_bac_taxonomy_file', help='file with bacterial GTDB taxonomy for genomes to consider (download from GTDB website)')
    parser.add_argument('gtdb_ar_taxonomy_file', help='file with archaeal GTDB taxonomy for genomes to consider (download from GTDB website)')
    parser.add_argument('silva_ssu_ref', help='FASTA file with SILVA Ref 16S sequences and taxonomy assignments')
    parser.add_argument('silva_lsu_ref', help='FASTA file with SILVA Ref 23S sequences and taxonomy assignments')
    parser.add_argument('ssu_blast_table', help='table with BLAST results between GTDB and SILVA Ref 16S rRNA sequences' )
    parser.add_argument('lsu_blast_table', help='table with BLAST results between GTDB and SILVA Ref 23S rRNA sequences' )
    parser.add_argument('output_dir', help="output directory")

    args = parser.parse_args()

    try:
        p = Create()
        p.run(args.gtdb_bac_taxonomy_file, 
                args.gtdb_ar_taxonomy_file,
                args.silva_ssu_ref,
                args.silva_lsu_ref,
                args.ssu_blast_table, 
                args.lsu_blast_table, 
                args.output_dir)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
