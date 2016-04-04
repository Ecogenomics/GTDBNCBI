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
    self.ssu_gg_taxonomy_file = os.path.join('ssu_gg_2013_08', 'ssu.taxonomy.tsv')
    self.ssu_silva_taxonomy_file = os.path.join('ssu_silva_199_gg_taxonomy', 'ssu.taxonomy.tsv')

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

  def _parse_taxonomy_file(self, genome_id, metadata_taxonomy_file, fout, prefix):
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
      header_line = f.readline() # consume header line
      if prefix not in self.taxonomy_headers:
        self.taxonomy_headers.add(prefix)

        fout.write('genome_id')
        headers = [prefix + '_' + x.strip().replace('ssu_', '') for x in header_line.split('\t')]
        fout.write('\t' + '\t'.join(headers) + '\n')

      # Report hit to longest 16S rRNA gene. It is possible that
      # the HMMs identified a putative 16S rRNA gene, but that
      # there was no valid BLAST hit.
      longest_query_len = 0
      longest_ssu_hit_info = None
      identified_ssu_genes = 0
      for line in f:
        identified_ssu_genes += 1
        line_split = line.strip().split('\t')
        query_len = int(line_split[2])
        if query_len > longest_query_len:
          longest_query_len = query_len
          longest_ssu_hit_info = line_split

      if longest_ssu_hit_info:
        fout.write(genome_id)
        fout.write('\t' + '\t'.join(longest_ssu_hit_info))
        fout.write('\n')

      return identified_ssu_genes

  def run(self, genbank_genome_dir, refseq_genome_dir, user_genome_dir, output_dir):
    """Create metadata tables."""

    fout_nt = open(os.path.join(output_dir, 'metadata_nt.tsv'), 'w')
    fout_gene = open(os.path.join(output_dir, 'metadata_gene.tsv'), 'w')
    fout_gg_taxonomy = open(os.path.join(output_dir, 'metadata_ssu.gg_2013_08.tsv'), 'w')
    fout_silva_taxonomy = open(os.path.join(output_dir, 'metadata_ssu.silva_199_gg_taxa.tsv'), 'w')
    fout_ssu_count = open(os.path.join(output_dir, 'metadata_ssu_count.tsv'), 'w')

    fout_ssu_count.write('%s\t%s\n' % ('genome_id', 'ssu_count'))

    # generate metadata for NCBI assemblies
    for ncbi_genome_dir in [genbank_genome_dir, refseq_genome_dir]:
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
            if len(processed_assemblies[accession]) >= 2:
              continue

            full_assembly_dir = os.path.join(full_species_dir, assembly_dir)

            protein_file = os.path.join(full_assembly_dir, assembly_dir + '_protein.faa')
            if not os.path.exists(protein_file):
              continue

            metadata_nt_file = os.path.join(full_assembly_dir, self.metadata_nt_file)
            self._parse_nt(accession, metadata_nt_file, fout_nt)

            metadata_gene_file = os.path.join(full_assembly_dir, self.metadata_gene_file)
            self._parse_gene(accession, metadata_gene_file, fout_gene)

            ssu_gg_taxonomy_file = os.path.join(full_assembly_dir, self.ssu_gg_taxonomy_file)
            self._parse_taxonomy_file(accession, ssu_gg_taxonomy_file, fout_gg_taxonomy, 'ssu_gg_2013_08')

            ssu_silva_taxonomy_file = os.path.join(full_assembly_dir, self.ssu_silva_taxonomy_file)
            ssu_count = self._parse_taxonomy_file(accession, ssu_silva_taxonomy_file, fout_silva_taxonomy, 'ssu_silva_199_gg_taxa')

            fout_ssu_count.write('%s\t%d\n' % (accession, ssu_count))

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
            self._parse_taxonomy_file(genome_id, ssu_gg_taxonomy_file, fout_gg_taxonomy, 'ssu_gg_2013_08')

            ssu_silva_taxonomy_file = os.path.join(full_genome_dir, self.ssu_silva_taxonomy_file)
            ssu_count = self._parse_taxonomy_file(genome_id, ssu_silva_taxonomy_file, fout_silva_taxonomy, 'ssu_silva_199_gg_taxa')

            fout_ssu_count.write('%s\t%d\n' % (genome_id, ssu_count))

    fout_nt.close()
    fout_gene.close()
    fout_gg_taxonomy.close()
    fout_silva_taxonomy.close()
    fout_ssu_count.close()

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
