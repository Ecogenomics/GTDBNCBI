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

__prog_name__ = 'silva_identify_bad_seqs.py'
__prog_desc__ = 'Parse tables with SILVA assignments to identify potentially erroneous 16S and 23S rRNA genes in GTDB genomes.'

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
    """Parse tables with SILVA assignments to identify potentially erroneous 16S and 23S rRNA genes in GTDB genomes."""

    def __init__(self):
        """Initialize."""
        
        pass
        
    def _canonical_taxonomy(self, silva_taxa):
        """Create canonical taxonomy with binomial species names."""
        
        if len(silva_taxa) >= 6:
            # assume last entry is a species designation
            sp = silva_taxa[-1].replace('Candidatus ', '').replace('uncultured ', '')
            sp_tokens = sp.split()
            
            if len(sp_tokens) > 2:
                sp = ' '.join(sp_tokens[0:2])
                
            silva_taxa[-1] = sp
        
        return silva_taxa
        
    def _multigene_silva_assignment_test(self, 
                                            rna_table, 
                                            gtdb_taxonomy,
                                            rna_ref,
                                            test_type, 
                                            test_note, 
                                            fout):
        """Check for congruent SILVA assignments for genomes with multiple 16S or 23S rRNA genes."""
        
        rna = {}
        with open(rna_table) as f:
            header = f.readline().strip().split('\t')
            
            silva_taxonomy_index = header.index('SILVA taxonomy')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gtdb_gene_id = line_split[0]
                gid = gtdb_gene_id.split('~')[0]
                
                silva_taxa = [t.strip() for t in line_split[silva_taxonomy_index].split(';')]
                canonical_silva_taxa = self._canonical_taxonomy(silva_taxa)
                
                if gid in rna:
                    prev_silva_taxa, prev_gtdb_gene_id = rna[gid]
                    
                    for rank_index, (prev_taxon, taxon) in enumerate(zip(prev_silva_taxa, canonical_silva_taxa)):
                        if prev_taxon != taxon:
                            fout.write('%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (gid,
                                                                                    test_type,
                                                                                    rank_index,
                                                                                    prev_taxon,
                                                                                    taxon,
                                                                                    '; '.join(prev_silva_taxa),
                                                                                    '; '.join(silva_taxa),
                                                                                    '; '.join(gtdb_taxonomy[gid]),
                                                                                    test_note,
                                                                                    prev_gtdb_gene_id in rna_ref,
                                                                                    gtdb_gene_id in rna_ref))
                                                                                    
                            if prev_gtdb_gene_id in rna_ref:
                                rna_len, contig_len = rna_ref[prev_gtdb_gene_id]
                                fout.write('\t%s\t%d\t%d' % (prev_gtdb_gene_id, rna_len, contig_len))
                            else:
                                fout.write('\t%s\t-\t-' % prev_gtdb_gene_id)
                                
                            if gtdb_gene_id in rna_ref:
                                rna_len, contig_len = rna_ref[gtdb_gene_id]
                                fout.write('\t%s\t%d\t%d' % (gtdb_gene_id, rna_len, contig_len))
                            else:
                                fout.write('\t%s\t-\t-' % gtdb_gene_id)
                                
                            fout.write('\n')
                            
                            break
                
                rna[gid] = (silva_taxa, gtdb_gene_id)
                
    def _read_rna_assignments(self, rna_table):
        """Read taxonomic assignments to rNA genes in a genome."""
        
        rna_assignments = defaultdict(lambda: {})
        with open(rna_table) as f:
            header = f.readline().strip().split('\t')
            
            silva_taxonomy_index = header.index('SILVA taxonomy')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gtdb_gene_id = line_split[0]
                gid = gtdb_gene_id.split('~')[0]
                
                silva_taxa = [t.strip() for t in line_split[silva_taxonomy_index].split(';')]
                
                rna_assignments[gid][gtdb_gene_id] = silva_taxa
                
        return rna_assignments
                
    def _ssu_lsu_silva_assignment_test(self, 
                                        ssu_table, 
                                        lsu_table, 
                                        gtdb_taxonomy, 
                                        ssu_ref,
                                        lsu_ref,
                                        test_type, 
                                        test_note, 
                                        fout):
        """Check is SSU and LSU genes in a genome have the same SILVA assignment."""
        
        ssu_assignments = self._read_rna_assignments(ssu_table)
        lsu_assignments = self._read_rna_assignments(lsu_table)
        
        for gid in ssu_assignments:
            if not gid in lsu_assignments:
                continue
                
            for ssu_gene_id in ssu_assignments[gid]:
                ssu_canonical_taxa = self._canonical_taxonomy(ssu_assignments[gid][ssu_gene_id])
                
                for lsu_gene_id in lsu_assignments[gid]:
                    lsu_canonical_taxa = self._canonical_taxonomy(lsu_assignments[gid][lsu_gene_id])
                
                    for rank_index, (ssu_taxon, lsu_taxon) in enumerate(zip(ssu_canonical_taxa, lsu_canonical_taxa)):
                            if ssu_taxon != lsu_taxon:
                                fout.write('%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (gid,
                                                                                test_type,
                                                                                rank_index,
                                                                                ssu_taxon,
                                                                                lsu_taxon,
                                                                                '; '.join(ssu_assignments[gid][ssu_gene_id]),
                                                                                '; '.join(lsu_assignments[gid][lsu_gene_id]),
                                                                                '; '.join(gtdb_taxonomy[gid]),
                                                                                test_note,
                                                                                ssu_gene_id in ssu_ref,
                                                                                lsu_gene_id in lsu_ref))
                                
                                if ssu_gene_id in ssu_ref:
                                    rna_len, contig_len = ssu_ref[ssu_gene_id]
                                    fout.write('\t%s\t%d\t%d' % (ssu_gene_id, rna_len, contig_len))
                                else:
                                    fout.write('\t%s\t-\t-' % ssu_gene_id)
                                    
                                if lsu_gene_id in lsu_ref:
                                    rna_len, contig_len = lsu_ref[lsu_gene_id]
                                    fout.write('\t%s\t%d\t%d' % (lsu_gene_id, rna_len, contig_len))
                                else:
                                    fout.write('\t%s\t-\t-' % lsu_gene_id)
                                fout.write('\n')
                                
                                break
                    
    def run(self, 
                gtdb_bac_taxonomy_file, 
                gtdb_ar_taxonomy_file, 
                ssu_silva_table, 
                lsu_silva_table,
                ssu_info_file,
                lsu_info_file,
                output_dir):
        """Parse tables with SILVA assignments to identify potentially erroneous 16S and 23S rRNA genes in GTDB genomes."""
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        fout = open(os.path.join(output_dir, 'silva_incongruence_test.tsv'), 'w')
        fout.write('Genome ID\tTest\tIncongruent rank')
        fout.write('\tSILVA taxon A\tSILVA taxon B\tSILVA taxonomy A\tSILVA taxonomy B')
        fout.write('\tGTDB taxonomy\tNote')
        fout.write('\tIn reference tree A\tIn reference tree B')
        fout.write('\tGene ID A\trRNA length\tContig length\tGene ID B\trRNA length\tContig length\n')
        
        # read GTDB taxonomy
        print('Reading GTDB taxonomy.')
        gtdb_bac_taxonomy = Taxonomy().read(gtdb_bac_taxonomy_file)
        gtdb_ar_taxonomy = Taxonomy().read(gtdb_ar_taxonomy_file)
        gtdb_taxonomy = gtdb_bac_taxonomy.copy()
        gtdb_taxonomy.update(gtdb_ar_taxonomy)
        
        # read genomes in SSU and LSU trees
        print('Reading genomes in 16S and 23S gene trees.')
        ssu_ref = {}
        with open(ssu_info_file) as f:
            header = f.readline().strip().split('\t')
            rna_length_index = header.index('SSU gene length')
            contig_len_index = header.index('Sequence length')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gene_id = line_split[0]
                contig_id = line_split[1]
                gene_id = gene_id.split('~')[0] + '~' + contig_id
                
                rna_length = int(line_split[rna_length_index])
                contig_length = int(line_split[contig_len_index])
                ssu_ref[gene_id] = (rna_length, contig_length)
                
        lsu_ref = {}
        with open(lsu_info_file) as f:
            header = f.readline().strip().split('\t')
            rna_length_index = header.index('SSU gene length')
            contig_len_index = header.index('Sequence length')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gene_id = line_split[0]
                contig_id = line_split[1]
                gene_id = gene_id.split('~')[0] + '~' + contig_id
                
                rna_length = int(line_split[rna_length_index])
                contig_length = int(line_split[contig_len_index])
                lsu_ref[gene_id] = (rna_length, contig_length)

        # run tests to find potentially incongruent 16S or 23S rRNA genes
        print('Performing tests to identify potentially incongruent 16S or 23S rRNA genes.')
        self._multigene_silva_assignment_test(ssu_silva_table, 
                                                gtdb_taxonomy, 
                                                ssu_ref,
                                                'SSU', 
                                                'Genome has multiple 16S rRNA genes with incongrent SILVA assignments.',
                                                fout)
                                        
        self._multigene_silva_assignment_test(lsu_silva_table, 
                                                gtdb_taxonomy, 
                                                lsu_ref,
                                                'LSU', 
                                                'Genome has multiple 23S rRNA genes with incongrent SILVA assignments.',
                                                fout)
                                        
        self._ssu_lsu_silva_assignment_test(ssu_silva_table,
                                            lsu_silva_table,
                                            gtdb_taxonomy,
                                            ssu_ref,
                                            lsu_ref,
                                            'SSU/LSU',
                                            'Genome has a 16S and 23S rRNA gene with incongruent SILVA assignments.',
                                            fout)
                                        
        fout.close()

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_bac_taxonomy_file', help='file with bacterial GTDB taxonomy for genomes to consider (download from GTDB website)')
    parser.add_argument('gtdb_ar_taxonomy_file', help='file with archaeal GTDB taxonomy for genomes to consider (download from GTDB website)')
    parser.add_argument('ssu_silva_table', help='file with 16S rRNA SILVA assignments for GTDB genomes')
    parser.add_argument('lsu_silva_table', help='file with 23S rRNA SILVA assignments for GTDB genomes')
    parser.add_argument('ssu_info_file', help='file with 16S rRNA information for genes in 16S rRNA gene tree')
    parser.add_argument('lsu_info_file', help='file with 23S rRNA information for genes in 23S rRNA gene tree')
    parser.add_argument('output_dir', help="output directory")

    args = parser.parse_args()

    try:
        p = Create()
        p.run(args.gtdb_bac_taxonomy_file, 
                args.gtdb_ar_taxonomy_file, 
                args.ssu_silva_table, 
                args.lsu_silva_table, 
                args.ssu_info_file,
                args.lsu_info_file,
                args.output_dir)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
