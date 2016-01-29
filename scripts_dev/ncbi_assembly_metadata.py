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

__prog_name__ = 'ncbi_assembly_metadata.py'
__prog_desc__ = 'Produce metadata file from NCBI assembly statistics.'

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
import argparse
from collections import defaultdict

from numpy import (zeros as np_zeros)

class GenericFeatureParser():
  """Parses generic feature file (GFF)."""
  def __init__(self, filename):
    self.cds_count = 0
    self.tRNA_count = 0
    self.rRNA_count = 0
    self.rRNA_16S_count = 0
    self.ncRNA_count = 0

    self.genes = {}
    self.last_coding_base = {}

    self._parse(filename)

    self.coding_mask = {}
    for seq_id in self.genes:
      self.coding_mask[seq_id] = self._coding_mask(seq_id)

  def _parse(self, gff_file):
    """Parse GFF file.

    Parameters
    ----------
    gff_file : str
      Generic feature file to parse.
    """

    for line in open(gff_file):
      if line[0] == '#':
        continue

      line_split = line.split('\t')
      if line_split[2] == 'tRNA':
        self.tRNA_count += 1
      elif line_split[2] == 'rRNA':
        self.rRNA_count += 1

        if 'product=16S ribosomal RNA' in line_split[8]:
          self.rRNA_16S_count += 1
      elif line_split[2] == 'ncRNA':
        self.ncRNA_count += 1
      elif line_split[2] == 'CDS':
        self.cds_count += 1

        seq_id = line_split[0]
        if seq_id not in self.genes:
          gene_count = 0
          self.genes[seq_id] = {}
          self.last_coding_base[seq_id] = 0

        gene_id = seq_id + '_' + str(gene_count)
        gene_count += 1

        start = int(line_split[3])
        end = int(line_split[4])

        self.genes[seq_id][gene_id] = [start, end]
        self.last_coding_base[seq_id] = max(self.last_coding_base[seq_id], end)

  def _coding_mask(self, seq_id):
    """Build mask indicating which bases in a sequences are coding."""

    # safe way to calculate coding bases as it accounts
    # for the potential of overlapping genes
    coding_mask = np_zeros(self.last_coding_base[seq_id])
    for pos in self.genes[seq_id].values():
      coding_mask[pos[0]:pos[1] + 1] = 1

    return coding_mask

  def coding_bases(self, seq_id):
    """Calculate number of coding bases in sequence."""

    # check if sequence has any genes
    if seq_id not in self.genes:
      return 0

    return np_sum(self.coding_mask[seq_id])

  def total_coding_bases(self):
    """Calculate total number of coding bases.

    Returns
    -------
    int
      Number of coding bases.
    """

    coding_bases = 0
    for seq_id in self.genes:
      coding_bases += self.coding_bases(seq_id)

    return int(coding_bases)

class Metadata(object):
  """Create metadata file from the assembly stats file of each NCBI assembly."""

  def __init__(self):
    self.fields = ['Assembly Name', 'Organism name', 'Taxid', 'Submitter', 'Date']
    self.fields.extend(['BioSample', 'Assembly type', 'Release type', 'Assembly level'])
    self.fields.extend(['Genome representation', 'GenBank Assembly Accession', 'RefSeq Assembly and GenBank Assemblies Identical'])

    self.stats = ['molecule-count', 'contig-count', 'scaffold-count', 'region-count', 'top-level-count']
    self.stats.extend(['spanned-gaps', 'total-gap-length', 'total-length', 'ungapped-length', 'unspanned-gaps'])
    self.stats.extend(['contig-N50', 'scaffold-N50', 'scaffold-L50', 'scaffold-N75', 'scaffold-N90'])

    self.gff_fields = ['protein_count', 'tRNA_count', 'ncRNA_count', 'rRNA_count', 'ssu_count']

    self.gbff_fields = ['translation_table']

    self.stats_info = {}

  def _parse_assembly_stats(self, assembly_stat_file):
    """Parse data from assembly stats file.

    Parameters
    ----------
    assembly_stat_file : str
      NCBI assembly stat file to parse.

    Returns
    -------
    list
      Parsed metadata in canonical order.
    """

    metadata_fields = [''] * len(self.fields)
    metadata_stats = [''] * len(self.stats)

    file_section = 'Assembly info'
    for line in open(assembly_stat_file):
      if 'Assembly Statistics Report' in line:
        file_section = 'ASR'
      elif 'Statistic Types' in line:
        file_section = 'ST'
      elif 'Sequence-type Description' in line:
        file_section = 'SD'

      if file_section == 'ASR' and ':' in line:
        field = line[2:line.find(':')]
        value = line[line.find(':')+1:].strip()

        if field in self.fields:
          if field == 'Organism name' and '(' in value:
            metadata_index = self.fields.index(field)
            metadata_fields[metadata_index] = value[0:value.find('(')].strip()
          else:
            metadata_index = self.fields.index(field)
            metadata_fields[metadata_index] = value
      elif file_section == 'ST':
        line_split = line.split('\t')
        if len(line_split) == 2:
          field = line_split[0][2:]
          desc = line_split[1].strip()
          if field in self.stats:
            self.stats_info[field] = desc
      elif file_section == 'SD':
        line_split = line.split()
        if len(line_split) == 6 and line_split[0] == 'all':
          field = line_split[4]
          value = line_split[5].strip()

          metadata_index = self.stats.index(field)
          metadata_stats[metadata_index] = value

    return metadata_fields, metadata_stats

  def _parse_gff(self, gff_file):
    """Parse statistics from generic feature file (GFF)."""

    metadata_gff = [''] * len(self.gff_fields)

    if not os.path.exists(gff_file):
      return metadata_gff

    gff_parser = GenericFeatureParser(gff_file)
    metadata_gff[self.gff_fields.index('protein_count')] = gff_parser.cds_count
    metadata_gff[self.gff_fields.index('tRNA_count')] = gff_parser.tRNA_count
    metadata_gff[self.gff_fields.index('ncRNA_count')] = gff_parser.ncRNA_count
    metadata_gff[self.gff_fields.index('rRNA_count')] = gff_parser.rRNA_count
    metadata_gff[self.gff_fields.index('ssu_count')] = gff_parser.rRNA_16S_count

    return metadata_gff

  def _parse_gbff(self, genbank_file):
    """Parse statistics from GenBank file."""

    metadata_gbff = [''] * len(self.gbff_fields)

    if not os.path.exists(genbank_file):
      return metadata_gbff

    for line in open(genbank_file):
      if '/transl_table=' in line:
        translation_table = line[line.rfind('=')+1:].strip()
        metadata_gbff[self.gbff_fields.index('translation_table')] = translation_table
        break

    return metadata_gbff

  def run(self, genome_dir, output_file):
    """Create metadata by parsing assembly stats files."""

    fout = open(output_file, 'w')
    fout.write('Assembly accession')
    fout.write('\t' + '\t'.join(['ncbi_' + x.lower().replace(' ', '_') for x in self.fields]))
    fout.write('\t' + '\t'.join(['ncbi_' + x.lower().replace('-', '_') for x in self.stats]))
    fout.write('\t' + '\t'.join(['ncbi_' + x.lower() for x in self.gff_fields]))
    fout.write('\t' + '\t'.join(['ncbi_' + x.lower() for x in self.gbff_fields]))
    fout.write('\n')

    processed_assemblies = defaultdict(list)
    for domain in ['archaea', 'bacteria']:
      domain_dir = os.path.join(genome_dir, domain)
      for species_dir in os.listdir(domain_dir):
        full_species_dir = os.path.join(domain_dir, species_dir)
        for assembly_dir in os.listdir(full_species_dir):
          accession = assembly_dir[0:assembly_dir.find('_', 4)]

          processed_assemblies[accession].append(species_dir)
          if len(processed_assemblies[accession]) >= 2:
            print '%s\t%s' % (accession, ','.join(processed_assemblies[accession]))
            continue

          full_assembly_dir = os.path.join(full_species_dir, assembly_dir)

          protein_file = os.path.join(full_assembly_dir, assembly_dir + '_protein.faa')
          if not os.path.exists(protein_file):
            continue

          assembly_stat_file = os.path.join(full_assembly_dir, assembly_dir + '_assembly_stats.txt')
          if os.path.exists(assembly_stat_file):
            metadata_fields, metadata_stats = self._parse_assembly_stats(assembly_stat_file)
            fout.write(accession + '\t%s\t%s' % ('\t'.join(metadata_fields),
                                                  '\t'.join(metadata_stats)))

            gff_file = os.path.join(full_assembly_dir, assembly_dir + '_genomic.gff')
            gff_stats = self._parse_gff(gff_file)
            fout.write('\t%s' % ('\t'.join(map(str, gff_stats))))

            genbank_file = os.path.join(full_assembly_dir, assembly_dir + '_genomic.gbff')
            gbff_stats = self._parse_gbff(genbank_file)
            fout.write('\t%s' % ('\t'.join(map(str, gbff_stats))))

          fout.write('\n')

    fout.close()

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('genome_dir', help='base directory leading to NCBI archaeal and bacterial genome assemblies')
  parser.add_argument('output_file', help='output metadata file')

  args = parser.parse_args()

  try:
    p = Metadata()
    p.run(args.genome_dir, args.output_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
