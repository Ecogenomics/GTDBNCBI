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
        self.lsu_silva_taxonomy_file = os.path.join('rna_silva', 'lsu_23S.taxonomy.tsv')
        self.lsu_silva_fna_file = os.path.join('rna_silva', 'lsu_23S.fna')
        self.lsu_silva_summary_file = os.path.join('rna_silva', 'lsu_23S.hmm_summary.tsv')

        self.write_nt_header = True
        self.write_gene_header = True
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
          Full path to file containing 16S rRNA metadata.
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

    def run(self, genbank_genome_dir, refseq_genome_dir, user_genome_dir, output_dir):
        """Create metadata tables."""

        fout_nt = open(os.path.join(output_dir, 'metadata_nt.tsv'), 'w')
        fout_gene = open(os.path.join(output_dir, 'metadata_gene.tsv'), 'w')
        fout_gg_taxonomy = open(os.path.join(output_dir, 'metadata_ssu_gg.tsv'), 'w')
        fout_ssu_silva_taxonomy = open(os.path.join(output_dir, 'metadata_ssu_silva.tsv'), 'w')
        fout_lsu_silva_taxonomy = open(os.path.join(output_dir, 'metadata_lsu_silva.tsv'), 'w')
        fout_ssu_silva_count = open(os.path.join(output_dir, 'metadata_ssu_silva_count.tsv'), 'w')
        fout_lsu_silva_count = open(os.path.join(output_dir, 'metadata_lsu_silva_count.tsv'), 'w')

        fout_ssu_silva_count.write('%s\t%s\n' % ('genome_id', 'ssu_count'))
        fout_lsu_silva_count.write('%s\t%s\n' % ('genome_id', 'lsu_count'))

        # generate metadata for NCBI assemblies
        # for ncbi_genome_dir in [genbank_genome_dir, refseq_genome_dir]:
        for ncbi_genome_dir in [genbank_genome_dir]:
            processed_assemblies = defaultdict(list)
            print 'Reading NCBI assembly directories: %s' % ncbi_genome_dir
            processed_assemblies = defaultdict(list)
            for domain in ['archaea', 'bacteria']:
                domain_dir = os.path.join(ncbi_genome_dir, domain)
                for species_dir in os.listdir(domain_dir):
                    full_species_dir = os.path.join(domain_dir, species_dir)
                    for assembly_dir in os.listdir(full_species_dir):
                        no_gc = ("GCA_000820785.2", "GCA_001059375.1", "GCA_001642695.1", "GCA_001642875.1", "GCA_001643385.1", "GCA_001643435.1", "GCA_001643455.1", "GCA_001643475.1", "GCA_001643485.1", "GCA_001643515.1", "GCA_001643535.1", "GCA_001643555.1", "GCA_001643565.1", "GCA_001643585.1", "GCA_001643615.1", "GCA_001643635.1", "GCA_001647515.1", "GCA_001649695.1", "GCA_001649735.1", "GCA_001649885.1", "GCA_001649915.1", "GCA_001655195.1", "GCA_001657295.1", "GCA_001657305.1", "GCA_001657315.1", "GCA_001657325.1", "GCA_001657375.1", "GCA_001657385.1", "GCA_001657395.1", "GCA_001657925.1", "GCA_001674955.1", "GCA_001678045.1", "GCA_001678065.1", "GCA_001683835.1", "GCA_001683845.1", "GCA_001683895.1", "GCA_001683905.1", "GCA_001683985.1", "GCA_001684155.1", "GCA_001684175.1", "GCA_001686345.1", "GCA_001686385.1", "GCA_001689405.1", "GCA_001689415.1", "GCA_001689425.1", "GCA_001689445.1", "GCA_001689485.1", "GCA_001689495.1", "GCA_001689515.1", "GCA_001689535.1", "GCA_001689565.1", "GCA_001689575.1", "GCA_001689585.1", "GCA_001689615.1", "GCA_001689645.1", "GCA_001689655.1", "GCA_001689665.1", "GCA_001689685.1", "GCA_001689725.1", "GCA_001700485.1", "GCA_001700545.1", "GCA_001701065.1", "GCA_001701075.1", "GCA_001701105.1", "GCA_001701115.1", "GCA_001701135.1", "GCA_001701165.1", "GCA_001701175.1", "GCA_001701195.1", "GCA_001701225.1", "GCA_001701235.1", "GCA_001701255.1", "GCA_001701285.1", "GCA_001701295.1", "GCA_001701305.1", "GCA_001702075.1", "GCA_001707145.1", "GCA_001707235.1", "GCA_001714685.1", "GCA_001717005.1", "GCA_001717015.1", "GCA_001717025.1", "GCA_001717035.1", "GCA_001717085.1", "GCA_001719265.1", "GCA_001719315.1", "GCA_001719375.1", "GCA_001719405.1", "GCA_001719445.1", "GCA_001719465.1", "GCA_001719545.1", "GCA_001723875.1", "GCA_001726005.1", "GCA_001726145.1", "GCA_001730645.1", "GCA_001735855.1", "GCA_001735875.1", "GCA_001735895.1", "GCA_001735915.1", "GCA_001742805.1", "GCA_001743105.1", "GCA_001743115.1", "GCA_001743125.1", "GCA_001743135.1", "GCA_001743185.1", "GCA_001743195.1", "GCA_001743215.1", "GCA_001743235.1", "GCA_001743265.1", "GCA_001743275.1", "GCA_001743285.1", "GCA_001743305.1", "GCA_001743345.1", "GCA_001743355.1", "GCA_001743375.1", "GCA_001743385.1", "GCA_001743425.1", "GCA_001743435.1", "GCA_001743455.1", "GCA_001743465.1", "GCA_001743495.1", "GCA_001746265.1", "GCA_001746285.1", "GCA_001746295.1", "GCA_001746305.1", "GCA_001746315.1", "GCA_001746365.1", "GCA_900002425.1", "GCA_900002435.1", "GCA_900002475.1", "GCA_900002485.1", "GCA_900002525.1", "GCA_900005695.1", "GCA_900006345.1", "GCA_900007725.1", "GCA_900007735.1", "GCA_900007745.1", "GCA_900007815.1", "GCA_900007825.1", "GCA_900007835.1", "GCA_900007845.1", "GCA_900007855.1", "GCA_900007865.1", "GCA_900007875.1", "GCA_900007885.1", "GCA_900007895.1", "GCA_900007905.1", "GCA_900007915.1", "GCA_900007925.1", "GCA_900007935.1", "GCA_900007945.1", "GCA_900007955.1", "GCA_900007965.1", "GCA_900007975.1", "GCA_900007985.1", "GCA_900007995.1", "GCA_900008005.1", "GCA_900008015.1", "GCA_900008025.1", "GCA_900008035.1", "GCA_900008045.1", "GCA_900008055.1", "GCA_900008065.1", "GCA_900008865.1", "GCA_900009155.1", "GCA_900009165.1", "GCA_900009175.1", "GCA_900009185.1", "GCA_900009195.1", "GCA_900009205.1", "GCA_900009285.1", "GCA_900009295.1", "GCA_900009305.1", "GCA_900009315.1", "GCA_900009565.1", "GCA_900009575.1", "GCA_900009585.1", "GCA_900009595.1", "GCA_900009605.1", "GCA_900009615.1", "GCA_900009625.1", "GCA_900009635.1", "GCA_900009645.1", "GCA_900009745.1", "GCA_900009755.1", "GCA_900009765.1", "GCA_900009775.1", "GCA_900009785.1", "GCA_900009835.1", "GCA_900009865.1", "GCA_900009875.1", "GCA_900009885.1", "GCA_900009895.1", "GCA_900009905.1", "GCA_900009915.1", "GCA_900009985.1", "GCA_900010075.1", "GCA_900010085.1", "GCA_900010095.1", "GCA_900010195.1", "GCA_900010245.1", "GCA_900010255.1", "GCA_900010265.1", "GCA_900010275.1", "GCA_900010285.1", "GCA_900010435.1", "GCA_900010445.1", "GCA_900010455.1", "GCA_900010555.1", "GCA_900010565.1", "GCA_900010735.1", "GCA_900010745.1", "GCA_900011035.1", "GCA_900011045.1", "GCA_900011205.1", "GCA_900011215.1", "GCA_900011445.1", "GCA_900011455.1", "GCA_900011465.1", "GCA_900011475.1", "GCA_900011485.1", "GCA_900011545.1", "GCA_900011555.1", "GCA_900011765.1", "GCA_900011775.1", "GCA_900012315.1", "GCA_900012395.1", "GCA_900012665.1", "GCA_900013145.1", "GCA_900013295.1", "GCA_900013305.1", "GCA_900013315.1", "GCA_900013565.1", "GCA_900013575.1", "GCA_900014835.1", "GCA_900015085.1", "GCA_900015985.1", "GCA_900015995.1", "GCA_900016005.1", "GCA_900016015.1", "GCA_900016025.1", "GCA_900016035.1", "GCA_900016045.1", "GCA_900016055.1", "GCA_900016065.1", "GCA_900016075.1", "GCA_900016085.1", "GCA_900016095.1", "GCA_900016105.1", "GCA_900016115.1", "GCA_900016125.1", "GCA_900016135.1", "GCA_900016365.1", "GCA_900016375.1", "GCA_900019265.1", "GCA_900064405.1", "GCA_900064415.1", "GCA_900064425.1", "GCA_900064775.1", "GCA_900073015.1", "GCA_900074625.1", "GCA_900074875.1", "GCA_900080205.1", "GCA_900086605.1", "GCA_900086705.1", "GCA_900087555.1", "GCA_900087565.1", "GCA_900087675.1", "GCA_900087705.1", "GCA_900087745.1", "GCA_900087775.1", "GCA_900088145.1", "GCA_900089525.1", "GCA_900089535.1", "GCA_900089545.1", "GCA_900089555.1", "GCA_900089565.1", "GCA_900089575.1", "GCA_900089785.1", "GCA_900089995.1", "GCA_900091325.1", "GCA_900092645.1", "GCA_900092655.1", "GCA_900095705.1", "GCA_900095825.1", "GCA_900095835.1", "GCA_900095845.1", "GCA_900095855.1", "GCA_900095875.1")

                        accession = assembly_dir[0:assembly_dir.find('_', 4)]

                        processed_assemblies[accession].append(species_dir)
                        if len(processed_assemblies[accession]) >= 2 and assembly_dir.startswith(no_gc):
                            print assembly_dir
                            print "processed assemblies"
                            continue

                        full_assembly_dir = os.path.join(full_species_dir, assembly_dir)

                        #protein_file = os.path.join(full_assembly_dir, assembly_dir + '_protein.faa')
                        genome_id = assembly_dir[0:assembly_dir.find('_', 4)]
                        protein_file = os.path.join(full_assembly_dir, "prodigal", genome_id + "_protein.faa")
                        if not os.path.exists(protein_file) and assembly_dir.startswith(no_gc):
                            print assembly_dir
                            print "protein file"
                            continue

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

                        lsu_silva_taxonomy_file = os.path.join(full_assembly_dir, self.lsu_silva_taxonomy_file)
                        lsu_silva_fna_file = os.path.join(full_assembly_dir, self.lsu_silva_fna_file)
                        lsu_silva_summary_file = os.path.join(full_assembly_dir, self.lsu_silva_summary_file)
                        lsu_count = self._parse_taxonomy_file(accession, lsu_silva_taxonomy_file, fout_lsu_silva_taxonomy, 'lsu_silva', lsu_silva_fna_file, lsu_silva_summary_file)

                        fout_ssu_silva_count.write('%s\t%d\n' % (accession, ssu_count))
                        fout_lsu_silva_count.write('%s\t%d\n' % (accession, lsu_count))

        # generate metadata for user genomes
        print 'Reading user genome directories.'
        if user_genome_dir != 'NONE':
            for user_id in os.listdir(user_genome_dir):
                full_user_dir = os.path.join(user_genome_dir, user_id)
                if not os.path.isdir(full_user_dir):
                    continue

                for genome_id in os.listdir(full_user_dir):
                    # for genome_id in ["U_65402", "U_65403", "U_65404", "U_65405", "U_65406", "U_65407", "U_65408", "U_65409", "U_65410"]:
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

                    lsu_silva_taxonomy_file = os.path.join(full_genome_dir, self.lsu_silva_taxonomy_file)
                    lsu_silva_fna_file = os.path.join(full_genome_dir, self.lsu_silva_fna_file)
                    lsu_silva_summary_file = os.path.join(full_genome_dir, self.lsu_silva_summary_file)
                    lsu_count = self._parse_taxonomy_file(genome_id, lsu_silva_taxonomy_file, fout_lsu_silva_taxonomy, 'lsu_silva', lsu_silva_fna_file, lsu_silva_summary_file)

                    fout_ssu_silva_count.write('%s\t%d\n' % (genome_id, ssu_count))
                    fout_lsu_silva_count.write('%s\t%d\n' % (genome_id, lsu_count))

        fout_nt.close()
        fout_gene.close()
        fout_gg_taxonomy.close()
        fout_ssu_silva_taxonomy.close()
        fout_lsu_silva_taxonomy.close()
        fout_ssu_silva_count.close()
        fout_lsu_silva_count.close()

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
        p.run(args.genbank_genome_dir, args.refseq_genome_dir, args.user_genome_dir, args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
