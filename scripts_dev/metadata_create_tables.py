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

__prog_name__ = 'metadata_create_tables.py'
__prog_desc__ = 'Create metadata tables for all NCBI and user genomes.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import ntpath
import argparse
from collections import defaultdict
from biolib.seq_io import read_fasta


class MetadataTable(object):
    """Create metadata table for all NCBI and user genomes.

    This script assumes the scripts metadata_generate.py
    and ssu.py have been run in order to create the
    required metadata. Four tables are generated which
    specific nucleotide derived, gene derived, and SSU
    derived metadata.
    """

    def __init__(self):
        self.metadata_nt_file = 'metadata.genome_nt.tsv'
        self.metadata_gene_file = 'metadata.genome_gene.tsv'
        self.ssu_gg_taxonomy_file = os.path.join('ssu_gg', 'ssu.taxonomy.tsv')
        self.ssu_gg_fna_file = os.path.join('ssu_gg', 'ssu.fna')
        self.ssu_silva_taxonomy_file = os.path.join('rna_silva', 'ssu.taxonomy.tsv')
        self.ssu_silva_fna_file = os.path.join('rna_silva', 'ssu.fna')
        self.ssu_silva_summary_file = os.path.join('rna_silva', 'ssu.hmm_summary.tsv')
        self.lsu_silva_23s_taxonomy_file = os.path.join('rna_silva', 'lsu_23S.taxonomy.tsv')
        self.lsu_silva_23s_fna_file = os.path.join('rna_silva', 'lsu_23S.fna')
        self.lsu_silva_23s_summary_file = os.path.join('rna_silva', 'lsu_23S.hmm_summary.tsv')
        
        self.lsu_5S_fna_file = os.path.join('lsu_5S', 'lsu_5S.fna')
        self.lsu_5S_summary_file = os.path.join('lsu_5S', 'lsu_5S.hmm_summary.tsv')

        self.write_nt_header = True
        self.write_gene_header = True
        self.write_trna_header = True
        self.write_lsu_5S_header = True
        self.taxonomy_headers = set()

    def _parse_nt(self, genome_id, metadata_nt_file, fout):
        """Parse metadata file with information derived from nucleotide sequences."""

        if not os.path.exists(metadata_nt_file):
            return

        if self.write_nt_header:
            self.write_nt_header = False

            fout.write('genome_id')
            for line in open(metadata_nt_file):
                line_split = line.split('\t')
                fout.write('\t' + line_split[0].strip())
            fout.write('\n')

        fout.write(genome_id)
        for line in open(metadata_nt_file):
            line_split = line.split('\t')
            fout.write('\t' + line_split[1].strip())
        fout.write('\n')

    def _parse_gene(self, genome_id, metadata_gene_file, fout):
        """Parse metadata file with information derived from called genes."""

        if not os.path.exists(metadata_gene_file):
            return

        if self.write_gene_header:
            self.write_gene_header = False

            fout.write('genome_id')
            for line in open(metadata_gene_file):
                line_split = line.split('\t')
                fout.write('\t' + line_split[0].strip())
            fout.write('\n')

        fout.write(genome_id)
        for line in open(metadata_gene_file):
            line_split = line.split('\t')
            fout.write('\t' + line_split[1].strip())
        fout.write('\n')

    def _parse_taxonomy_file(self, genome_id, metadata_taxonomy_file, fout, prefix, fna_file, summary_file=None):
        """Parse metadata file with taxonomic information for 16S rRNA genes.

        Parameters
        ----------
        genome_id : str
          Unique identifier of genome.
        metadata_taxonomy_file : str
          Full path to file containing rRNA metadata.
        fout : file
          Output stream to populate with metadata.
        Prefix : str
          Prefix to append to metadata fields.

        Returns
        -------
        int
          Number of 16S rRNA genes identified in genome.
        """

        if not os.path.exists(metadata_taxonomy_file):
            return 0

        with open(metadata_taxonomy_file) as f:
            header_line = f.readline()  # consume header line
            if prefix not in self.taxonomy_headers:
                self.taxonomy_headers.add(prefix)

                fout.write('genome_id')
                headers = [prefix + '_' + x.strip().replace('ssu_', '') for x in header_line.split('\t')]
                fout.write('\t' + '\t'.join(headers))
                fout.write('\t{0}_sequence\t{0}_contig_len\n'.format(prefix))

            # Check the CheckM headers are consistent
            split_headers = header_line.rstrip().split("\t")
            for pos in range(0, len(split_headers)):
                header = split_headers[pos]
                if header == 'query_id':
                    query_id_pos = pos
                    break

            # Report hit to longest 16S rRNA gene. It is possible that
            # the HMMs identified a putative 16S rRNA gene, but that
            # there was no valid BLAST hit.
            longest_query_len = 0
            longest_ssu_hit_info = None
            identified_ssu_genes = 0
            for line in f:
                line_split = line.strip().split('\t')
                query_len = int(line_split[2])
                if query_len > longest_query_len:
                    longest_query_len = query_len
                    longest_ssu_hit_info = line_split
                    ssu_query_id = line_split[query_id_pos]

            if longest_ssu_hit_info:
                fout.write(genome_id)
                fout.write('\t' + '\t'.join(longest_ssu_hit_info))
                all_genes_dict = read_fasta(fna_file, False)
                sequence = all_genes_dict[ssu_query_id]
                fout.write('\t{0}'.format(sequence))
                if summary_file is not None and os.path.exists(summary_file):
                    with open(summary_file) as fsum:
                        header_line = fsum.readline()  # consume header line
                        header_list = [x.strip() for x in header_line.split('\t')]
                        idx_seq = header_list.index("Sequence length")
                        for line in fsum:
                            identified_ssu_genes += 1
                            sum_list = [x.strip() for x in line.split('\t')]
                            if sum_list[0] == ssu_query_id:
                                fout.write("\t{0}".format(sum_list[idx_seq]))

                fout.write('\n')

            return identified_ssu_genes
            
    def _parse_lsu_5S_files(self, accession, fout, fna_file, summary_file):
        """Parse information from 5S LSU files."""
        
        # check if a 5S sequence was identified
        if not os.path.exists(fna_file):
            return 0
        
        # write header
        if self.write_lsu_5S_header:
            fout.write('genome_id\tlsu_5s_query_id\tlsu_5s_length\tlsu_5s_contig_len\tlsu_5s_sequence\n')
            self.write_lsu_5S_header = False

        seqs = read_fasta(fna_file)
            
        identified_genes = 0
        longest_seq = 0
        longest_seq_id = None
        longest_contig_len = None
        if os.path.exists(summary_file):
            with open(summary_file) as fsum:
                header_line = fsum.readline()  # consume header line
                header_list = [x.strip() for x in header_line.split('\t')]
                idx_seq_len = header_list.index("Sequence length")
                for line in fsum:
                    identified_genes += 1
                    
                    line_split = map(str.strip, line.strip().split('\t'))
                    seq_id = line_split[0]
                    contig_len = int(line_split[idx_seq_len])
                    seq_len = len(seqs[seq_id])
                    
                    if seq_len > longest_seq:
                        longest_seq_id = seq_id
                        longest_seq = seq_len
                        longest_contig_len = contig_len
                        
        if longest_seq_id:
            fout.write('%s\t%s\t%d\t%d\t%s\n' % (accession, longest_seq_id, longest_seq, longest_contig_len, seqs[longest_seq_id]))
                        
        return identified_genes
            
    def _parse_trna_file(self, genome_id, trna_file, fout_trna_count):
        """Parse tRNA information."""
        
        if not os.path.exists(trna_file):
            return
        
        # write header
        if self.write_trna_header:
            fout_trna_count.write('genome_id\ttrna_count\ttrna_aa_count\ttrna_selenocysteine_count\n')
            self.write_trna_header = False
            
        # parse tRNA summary file
        trna_count = 0
        trna_selenocysteine_count = 0
        trna_aa_count = 0
        read_aa = False
        for line in open(trna_file):
            if line.startswith('tRNAs decoding Standard 20 AA'):
                trna_count = int(line.split(':')[1])
            elif line.startswith('Selenocysteine tRNAs (TCA)'):
                trna_selenocysteine_count = int(line.split(':')[1])
                trna_count += trna_selenocysteine_count
            elif line.startswith('Isotype / Anticodon Counts:'):
                read_aa = True
                
            if read_aa:
                # Parsing lines with the following format:
                # Ala   : 2	  AGC:         GGC:         CGC: 1       TGC: 1 
                if ' : ' in line:
                    line_split = line.split('\t')[0]
                    line_split = map(str.strip, line_split.split(':'))
                    if len(line_split[0]) == 3: # this is an amino acid
                        if int(line_split[1]) > 0:
                            trna_aa_count += 1

        fout_trna_count.write('%s\t%d\t%d\t%d\n' % (genome_id, trna_count, trna_aa_count, trna_selenocysteine_count))

    def run(self, genbank_genome_dir, refseq_genome_dir, user_genome_dir, output_dir):
        """Create metadata tables."""

        fout_nt = open(os.path.join(output_dir, 'metadata_nt.tsv'), 'w')
        fout_gene = open(os.path.join(output_dir, 'metadata_gene.tsv'), 'w')
        fout_gg_taxonomy = open(os.path.join(output_dir, 'metadata_ssu_gg.tsv'), 'w')
        fout_ssu_silva_taxonomy = open(os.path.join(output_dir, 'metadata_ssu_silva.tsv'), 'w')
        fout_lsu_silva_23s_taxonomy = open(os.path.join(output_dir, 'metadata_lsu_silva_23s.tsv'), 'w')
        fout_lsu_5S = open(os.path.join(output_dir, 'metadata_lsu_5S.tsv'), 'w')
        fout_ssu_silva_count = open(os.path.join(output_dir, 'metadata_ssu_silva_count.tsv'), 'w')
        fout_lsu_silva_23s_count = open(os.path.join(output_dir, 'metadata_lsu_silva_23s_count.tsv'), 'w')
        fout_lsu_5S_count = open(os.path.join(output_dir, 'metadata_lsu_5S_count.tsv'), 'w')
        fout_trna_count = open(os.path.join(output_dir, 'metadata_trna_count.tsv'), 'w')

        fout_ssu_silva_count.write('%s\t%s\n' % ('genome_id', 'ssu_count'))
        fout_lsu_silva_23s_count.write('%s\t%s\n' % ('genome_id', 'lsu_23s_count'))
        fout_lsu_5S_count.write('%s\t%s\n' % ('genome_id', 'lsu_5S_count'))

        # generate metadata for NCBI assemblies
        # for ncbi_genome_dir in [genbank_genome_dir, refseq_genome_dir]:
        for ncbi_genome_dir in [genbank_genome_dir, refseq_genome_dir]:
            if ncbi_genome_dir == 'NONE':
                continue
                
            processed_assemblies = defaultdict(list)
            print 'Reading NCBI assembly directories: %s' % ncbi_genome_dir
            processed_assemblies = defaultdict(list)
            for domain in ['archaea', 'bacteria']:
                domain_dir = os.path.join(ncbi_genome_dir, domain)
                for species_dir in os.listdir(domain_dir):
                    full_species_dir = os.path.join(domain_dir, species_dir)
                    for assembly_dir in os.listdir(full_species_dir):
                        accession = assembly_dir[0:assembly_dir.find('_', 4)]
                        processed_assemblies[accession].append(species_dir)
                        full_assembly_dir = os.path.join(full_species_dir, assembly_dir)

                        #protein_file = os.path.join(full_assembly_dir, assembly_dir + '_protein.faa')
                        genome_id = assembly_dir[0:assembly_dir.find('_', 4)]
                        protein_file = os.path.join(full_assembly_dir, "prodigal", genome_id + "_protein.faa")

                        metadata_nt_file = os.path.join(full_assembly_dir, self.metadata_nt_file)
                        self._parse_nt(accession, metadata_nt_file, fout_nt)

                        metadata_gene_file = os.path.join(full_assembly_dir, self.metadata_gene_file)
                        self._parse_gene(accession, metadata_gene_file, fout_gene)

                        ssu_gg_taxonomy_file = os.path.join(full_assembly_dir, self.ssu_gg_taxonomy_file)
                        ssu_gg_fna_file = os.path.join(full_assembly_dir, self.ssu_gg_fna_file)
                        self._parse_taxonomy_file(accession, ssu_gg_taxonomy_file, fout_gg_taxonomy, 'ssu_gg', ssu_gg_fna_file)

                        ssu_silva_taxonomy_file = os.path.join(full_assembly_dir, self.ssu_silva_taxonomy_file)
                        ssu_silva_fna_file = os.path.join(full_assembly_dir, self.ssu_silva_fna_file)
                        ssu_silva_summary_file = os.path.join(full_assembly_dir, self.ssu_silva_summary_file)
                        ssu_count = self._parse_taxonomy_file(accession, ssu_silva_taxonomy_file, fout_ssu_silva_taxonomy, 'ssu_silva', ssu_silva_fna_file, ssu_silva_summary_file)

                        lsu_silva_23s_taxonomy_file = os.path.join(full_assembly_dir, self.lsu_silva_23s_taxonomy_file)
                        lsu_silva_23s_fna_file = os.path.join(full_assembly_dir, self.lsu_silva_23s_fna_file)
                        lsu_silva_23s_summary_file = os.path.join(full_assembly_dir, self.lsu_silva_23s_summary_file)
                        lsu_23s_count = self._parse_taxonomy_file(accession, lsu_silva_23s_taxonomy_file, fout_lsu_silva_23s_taxonomy, 'lsu_silva_23s', lsu_silva_23s_fna_file, lsu_silva_23s_summary_file)

                        lsu_5S_fna_file = os.path.join(full_assembly_dir, self.lsu_5S_fna_file)
                        lsu_5S_summary_file = os.path.join(full_assembly_dir, self.lsu_5S_summary_file)
                        lsu_5S_count = self._parse_lsu_5S_files(accession, fout_lsu_5S, lsu_5S_fna_file, lsu_5S_summary_file)
                        
                        fout_ssu_silva_count.write('%s\t%d\n' % (accession, ssu_count))
                        fout_lsu_silva_23s_count.write('%s\t%d\n' % (accession, lsu_23s_count))
                        fout_lsu_5S_count.write('%s\t%d\n' % (accession, lsu_5S_count))
                        
                        trna_file = os.path.join(full_assembly_dir, 'trna', accession + '_trna_stats.tsv')
                        self._parse_trna_file(accession, trna_file, fout_trna_count)

        # generate metadata for user genomes
        print 'Reading user genome directories.'
        if user_genome_dir != 'NONE':
            for user_id in os.listdir(user_genome_dir):
                full_user_dir = os.path.join(user_genome_dir, user_id)
                if not os.path.isdir(full_user_dir):
                    continue

                for genome_id in os.listdir(full_user_dir):
                    full_genome_dir = os.path.join(full_user_dir, genome_id)

                    metadata_nt_file = os.path.join(full_genome_dir, self.metadata_nt_file)
                    self._parse_nt(genome_id, metadata_nt_file, fout_nt)

                    metadata_gene_file = os.path.join(full_genome_dir, self.metadata_gene_file)
                    self._parse_gene(genome_id, metadata_gene_file, fout_gene)

                    ssu_gg_taxonomy_file = os.path.join(full_genome_dir, self.ssu_gg_taxonomy_file)
                    ssu_gg_fna_file = os.path.join(full_genome_dir, self.ssu_gg_fna_file)
                    self._parse_taxonomy_file(genome_id, ssu_gg_taxonomy_file, fout_gg_taxonomy, 'ssu_gg', ssu_gg_fna_file)

                    ssu_silva_taxonomy_file = os.path.join(full_genome_dir, self.ssu_silva_taxonomy_file)
                    ssu_silva_fna_file = os.path.join(full_genome_dir, self.ssu_silva_fna_file)
                    ssu_silva_summary_file = os.path.join(full_genome_dir, self.ssu_silva_summary_file)
                    ssu_count = self._parse_taxonomy_file(genome_id, ssu_silva_taxonomy_file, fout_ssu_silva_taxonomy, 'ssu_silva', ssu_silva_fna_file, ssu_silva_summary_file)

                    lsu_silva_23s_taxonomy_file = os.path.join(full_genome_dir, self.lsu_silva_23s_taxonomy_file)
                    lsu_silva_23s_fna_file = os.path.join(full_genome_dir, self.lsu_silva_23s_fna_file)
                    lsu_silva_23s_summary_file = os.path.join(full_genome_dir, self.lsu_silva_23s_summary_file)
                    lsu_23s_count = self._parse_taxonomy_file(genome_id, lsu_silva_23s_taxonomy_file, fout_lsu_silva_23s_taxonomy, 'lsu_silva_23s', lsu_silva_23s_fna_file, lsu_silva_23s_summary_file)

                    lsu_5S_fna_file = os.path.join(full_genome_dir, self.lsu_5S_fna_file)
                    lsu_5S_summary_file = os.path.join(full_genome_dir, self.lsu_5S_summary_file)
                    lsu_5S_count = self._parse_lsu_5S_files(genome_id, fout_lsu_5S, lsu_5S_fna_file, lsu_5S_summary_file)
                    
                    fout_ssu_silva_count.write('%s\t%d\n' % (genome_id, ssu_count))
                    fout_lsu_silva_23s_count.write('%s\t%d\n' % (genome_id, lsu_23s_count))
                    fout_lsu_5S_count.write('%s\t%d\n' % (genome_id, lsu_5S_count))
                    
                    trna_file = os.path.join(full_genome_dir, 'trna', genome_id + '_trna_stats.tsv')
                    self._parse_trna_file(genome_id, trna_file, fout_trna_count)

        fout_nt.close()
        fout_gene.close()
        fout_gg_taxonomy.close()
        fout_ssu_silva_taxonomy.close()
        fout_lsu_silva_23s_taxonomy.close()
        fout_lsu_5S.close()
        fout_ssu_silva_count.close()
        fout_lsu_silva_23s_count.close()
        fout_lsu_5S_count.close()
        fout_trna_count.close()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genbank_genome_dir', help='base directory leading to NCBI GenBank archaeal and bacterial genome assemblies')
    parser.add_argument('refseq_genome_dir', help='base directory leading to NCBI RefSeq archaeal and bacterial genome assemblies')
    parser.add_argument('user_genome_dir', help='base directory leading to user genomes or NONE to skip')
    parser.add_argument('output_dir', help='output directory')

    args = parser.parse_args()

    try:
        p = MetadataTable()
        p.run(args.genbank_genome_dir, 
                args.refseq_genome_dir, 
                args.user_genome_dir, 
                args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
