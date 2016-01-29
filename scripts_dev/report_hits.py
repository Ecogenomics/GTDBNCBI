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

__prog_name__ = 'report_hits.py'
__prog_desc__ = 'Report hits to a specified gene across genomes.'

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


class ReportHits(object):
  """Report top hits to specific gene across all genomes."""

  def __init__(self):
    pass

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

  def run(self, ncbi_genome_dir, gene_id):
    """Report hits."""

    # parse top hit tables
    for domain in ['archaea', 'bacteria']:
      num_genomes = 0
      num_hits = 0

      domain_dir = os.path.join(ncbi_genome_dir, domain)
      for species_dir in os.listdir(domain_dir):
        full_species_dir = os.path.join(domain_dir, species_dir)
        for assembly_dir in os.listdir(full_species_dir):
          accession = assembly_dir[0:assembly_dir.find('_', 4)]

          full_assembly_dir = os.path.join(full_species_dir, assembly_dir)

          if 'PFAM' in gene_id:
            top_hit_file = os.path.join(full_assembly_dir, assembly_dir + '_pfam_tophit.tsv')
          else:
            top_hit_file = os.path.join(full_assembly_dir, assembly_dir + '_tigrfam_tophit.tsv')

          if os.path.exists(top_hit_file):
            num_genomes += 1
            for line in open(top_hit_file):
              line_split = line.split('\t')
              if gene_id in line_split[1]:
                #print accession, line_split[1].strip()
                num_hits += 1

      print 'Num. hits in %s: %d' % (domain, num_hits)
      print 'Num. genomes in %s: %d' % (domain, num_genomes)

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('ncbi_genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies')
  parser.add_argument('gene_id', help='Gene of interest (PFAM or TIGRFAM identifier)')

  args = parser.parse_args()

  try:
    p = ReportHits()
    p.run(args.ncbi_genome_dir, args.gene_id)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
