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

__prog_name__ = 'map_taxonomy_img_to_ncbi.py'
__prog_desc__ = 'Map genome tree taxonomy defined over IMG genomes to NCBI genomes.'

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

class MapTaxonomy(object):
  """Map taxonomic information specified for IMG genomes to their NCBI counter-parts.

  This task is made complicated as there is no easily obtained unique identifier for
  mapping NCBI genomes to IMG genomes. As such, the organism name in conjunction with
  the genome size is used to identify the same genomes between these two datasets.
  """

  def __init__(self):
    pass

  def run(self, ncbi_metadata_file, img_metadata_file, taxonomy_file, output_taxonomy_file):
    """Map taxonomic information through organism names and genome size."""

    # get map between organism names and NCBI assembly identifiers
    org_name_id_ncbi_assemblies = defaultdict(list)
    num_ncbi_genomes = 0
    with open(ncbi_metadata_file) as f:
      headers = [x.strip() for x in f.readline().split('\t')]
      organism_name_index = headers.index('ncbi_organism_name')
      genome_size_index = headers.index('ncbi_total_length')

      for line in f:
        line_split = [x.strip() for x in line.split('\t')]

        org_name_id = ' '.join(line_split[organism_name_index].split(' ')[0:2]) + '_' + line_split[genome_size_index]
        org_name_id_ncbi_assemblies[org_name_id].append(line_split[0])
        num_ncbi_genomes += 1

    # get map between IMG genome identifiers and organism name
    img_organism_name = {}
    with open(img_metadata_file) as f:
      headers = [x.strip() for x in f.readline().split('\t')]
      organism_name_index = headers.index('Genome Name / Sample Name')
      genome_size_index = headers.index('Genome Size')

      for line in f:
        line_split = [x.strip() for x in line.split('\t')]

        org_name_id = ' '.join(line_split[organism_name_index].split(' ')[0:2]) + '_' + line_split[genome_size_index]
        img_organism_name[line_split[0]] = org_name_id

    print 'IMG genomes: %d' % len(img_organism_name)

    # get genome tree taxonomy string for IMG and user genomes
    taxonomy = {}
    img_genomes_with_taxonomy = 0
    user_genome_with_taxonomy = 0
    for line in open(taxonomy_file):
      line_split = [x.strip() for x in line.split('\t')]
      taxonomy[line_split[0].replace('IMG_', '')] = line_split[1]

      if line_split[0].startswith('IMG_'):
        img_genomes_with_taxonomy += 1
      else:
        user_genome_with_taxonomy += 1

    print 'IMG and user genomes with taxonomy: %d, %d' % (img_genomes_with_taxonomy, user_genome_with_taxonomy)

    # get 'best' map between organism names and genome tree taxonomy
    if False:
      organism_name_taxonomy = {}
      for genome_id, organism_name in img_organism_name.iteritems():
        if genome_id not in taxonomy:
          continue

        taxonomy_str = taxonomy[genome_id]

        if organism_name not in organism_name_taxonomy:
          organism_name_taxonomy[organism_name] = taxonomy_str
        elif taxonomy_str != organism_name_taxonomy[organism_name]:
          taxa = taxonomy_str.split(';')
          if taxa[5][3:] + ' ' + taxa[6][3:] == organism_name:
            # matches organism name so is likely a better taxonomy identifier
            print organism_name
            print organism_name_taxonomy[organism_name]
            print taxonomy_str
            organism_name_taxonomy[organism_name] = taxonomy_str
            raw_input('*')


    # create new taxonomy file for NCBI assemblies
    missing_img_metadata = 0
    num_ncbi_assigned_taxonomy = 0
    num_user_genomes_assigned_taxonomy = 0

    processed_org_name_id = {}

    fout = open(output_taxonomy_file, 'w')
    for genome_id, taxonomy_str in taxonomy.iteritems():
      if genome_id.startswith('U_'):
        fout.write('%s\t%s\n' % (genome_id, taxonomy_str))
        num_user_genomes_assigned_taxonomy += 1
      else:
        org_name_id = img_organism_name.get(genome_id, None)
        if org_name_id in processed_org_name_id:
          if processed_org_name_id[org_name_id] != taxonomy_str:
            print org_name_id
            print processed_org_name_id[org_name_id]
            print taxonomy_str
            raw_input('****')
          continue

        if not org_name_id:
          missing_img_metadata += 1
        else:
          processed_org_name_id[org_name_id] = taxonomy_str

          for ncbi_assembly_id in org_name_id_ncbi_assemblies.get(org_name_id, []):
            if ncbi_assembly_id.startswith('GCA_'):
              ncbi_assembly_id = 'GB_' + ncbi_assembly_id
            elif ncbi_assembly_id.startswith('GCF_'):
              ncbi_assembly_id = 'RS_' + ncbi_assembly_id

            fout.write('%s\t%s\n' % (ncbi_assembly_id, taxonomy_str))
            num_ncbi_assigned_taxonomy += 1

    fout.close()

    print ''
    print 'Missing IMG metadata: %d' % missing_img_metadata
    print ''
    print 'NCBI genomes: %d' % num_ncbi_genomes
    print 'NCBI genomes assigned taxonomy: %d' % num_ncbi_assigned_taxonomy
    print 'User genomes assigned taxonomy: %d' % num_user_genomes_assigned_taxonomy

if __name__ == '__main__':
  print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
  print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('ncbi_metadata_file', help='metadata for NCBI genomes including taxon Id')
  parser.add_argument('img_metadata_file', help='metadata for IMG genomes including taxon id')
  parser.add_argument('taxonomy_file', help='taxonomy file with taxonomic strings for IMG and user genomes')
  parser.add_argument('output_taxonomy_file', help='new taxonomy file with taxonomic strings for NCBI assemblies and user genomes')

  args = parser.parse_args()

  try:
    p = MapTaxonomy()
    p.run(args.ncbi_metadata_file, args.img_metadata_file, args.taxonomy_file, args.output_taxonomy_file)
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise
