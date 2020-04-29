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

__prog_name__ = 'ncbi_taxonomy.py'
__prog_desc__ = 'Parse NCBI taxonomy files to produce simplified summary files.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.5'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import traceback
from collections import namedtuple, defaultdict

from biolib.taxonomy import Taxonomy


class TaxonomyNCBI(object):
    """Parse NCBI taxonomy files to produce a simplified summary file."""

    def __init__(self):
        self.NodeRecord = namedtuple(
            'NodeRecord', 'parent_tax_id rank division_id genetic_code_id')
        self.NameRecord = namedtuple('NamesRecord', 'name_txt')

        self.bacterial_division = '0'
        self.unassigned_division = '8'

    def _assembly_organism_name(self, refseq_archaea_assembly_file, refseq_bacteria_assembly_file,
                                genbank_archaea_assembly_file, genbank_bacteria_assembly_file, output_organism_name_file):
        """Parse out organism name for each genome."""

        fout = open(output_organism_name_file, 'w')
        for assembly_file in [refseq_archaea_assembly_file, refseq_bacteria_assembly_file,
                              genbank_archaea_assembly_file, genbank_bacteria_assembly_file]:
            with open(assembly_file) as f:
                f.readline()
                header = f.readline().strip().split('\t')
                org_name_index = header.index('organism_name')

                for line in f:
                    line_split = line.strip().split('\t')

                    gid = line_split[0]
                    if gid.startswith('GCA_'):
                        gid = 'GB_' + gid
                    else:
                        gid = 'RS_' + gid
                    org_name = line_split[org_name_index]
                    fout.write('%s\t%s\n' % (gid, org_name))
        fout.close()

    def _assembly_to_tax_id(self, refseq_archaea_assembly_file, refseq_bacteria_assembly_file,
                            genbank_archaea_assembly_file, genbank_bacteria_assembly_file):
        """Determine taxonomic identifier for each assembly.

        Returns
        -------
        dict : d[assembly_accession] -> tax_id
          Taxonomic identifier for each assembly.
        """

        d = {}
        for assembly_file in [refseq_archaea_assembly_file, refseq_bacteria_assembly_file,
                              genbank_archaea_assembly_file, genbank_bacteria_assembly_file, ]:
            with open(assembly_file) as f:
                headers = f.readline().strip().split('\t')
                try:
                    taxid_index = headers.index('taxid')
                except:
                    # look for taxid on the next line as NCBI sometimes puts
                    # an extra comment on the first line
                    headers = f.readline().split('\t')
                    taxid_index = headers.index('taxid')

                for line in f:
                    line_split = line.strip().split('\t')
                    assembly_accession = line_split[0]
                    taxid = line_split[taxid_index]

                    if assembly_accession in d:
                        print '[Error] Duplicate assembly accession: %s' % assembly_accession
                        sys.exit(-1)

                    d[assembly_accession] = taxid

        return d

    def _read_nodes(self, nodes_file):
        """Read NCBI nodes.dmp file.

        Parameters
        ----------
        nodes_file : str
          Path to NCBI nodes.dmp file.

        Returns
        -------
        dict : d[tax_id] -> NodeRecord
          Node record for all nodes.
        """

        d = {}
        for line in open(nodes_file):
            line_split = [t.strip() for t in line.split('|')]

            tax_id = line_split[0]
            parent_tax_id = line_split[1]
            rank = line_split[2]
            division_id = line_split[4]
            genetic_code_id = line_split[6]

            d[tax_id] = self.NodeRecord(
                parent_tax_id, rank, division_id, genetic_code_id)

        return d

    def _read_names(self, names_file):
        """Read NCBI names.dmp file.

        Parameters
        ----------
        names_file : str
          Path to NCBI names.dmp file.

        Returns
        ------- 
        dict : d[tax_id] -> NameRecord
          Name record of nodes marked as 'scientific name'.
        """

        d = {}
        for line in open(names_file):
            line_split = [t.strip() for t in line.split('|')]

            tax_id = line_split[0]
            name_txt = line_split[1]
            unique_name = line_split[2]
            name_class = line_split[3]

            if name_class == 'scientific name':
                d[tax_id] = self.NameRecord(name_txt)

        return d

    def _valid_species_name(self, species_name, require_full=True, require_prefix=True):
        """Check if species name is a valid binomial name."""

        if species_name == 's__':
            return True, None

        # remove single quotes as sometimes given for
        # candidatus species names
        species_name = species_name.replace("'", "")

        # test for prefix
        if require_prefix:
            if not species_name.startswith('s__'):
                return False, 'name is missing the species prefix'

        # remove prefix before testing other properties
        test_name = species_name
        if test_name.startswith('s__'):
            test_name = test_name[3:]

        # test for full name
        if require_full:
            if 'candidatus' in test_name.lower():
                if len(test_name.split(' ')) <= 2:
                    return False, 'name appears to be missing the generic name'
            else:
                if len(test_name.split(' ')) <= 1:
                    return False, 'name appears to be missing the generic name'

        # get putative binomial name
        if 'candidatus' in test_name.lower():
            sp_name = ' '.join(test_name.split()[0:3])
        else:
            sp_name = ' '.join(test_name.split()[0:2])

        # check for tell-tale signs on invalid species names
        if sp_name[0].islower():
            return False, 'first letter of name is lowercase'
        if sp_name.split()[-1].isupper():
            return False, 'first letter of specific name is uppercase'
        if " bacterium" in sp_name.lower():
            return False, "name contains the word 'bacterium'"
        if " bacteirum" in sp_name.lower():
            return False, "name contains the word 'bacteirum'"
        if " bacteria" in sp_name.lower():
            return False, "name contains the word 'bacteria'"
        if " archaea" in sp_name.lower():
            return False, "name contains the word 'archaea'"
        if " archaeon" in sp_name.lower():
            return False, "name contains the word 'archaeon'"
        if " archeaon" in sp_name.lower():
            return False, "name contains the word 'archeaon'"
        if " archaeum" in sp_name.lower():
            return False, "name contains the word 'archaeum'"
        if "cyanobacterium" in sp_name.lower().split()[-1]:
            return False, "specific name is 'cyanobacterium'"
        if " group" in sp_name.lower():
            return False, "name contains 'group'"
        if " subdivision" in sp_name.lower():
            return False, "name contains 'subdivision'"
        if " taxon" in sp_name.lower():
            return False, "name contains 'taxon'"
        if " cluster" in sp_name.lower():
            return False, "name contains 'cluster'"
        if " clade" in sp_name.lower():
            return False, "name contains 'clade'"
        if " of " in sp_name.lower():
            return False, "name contains 'of'"
        if 'sp.' in sp_name.lower():
            return False, "name contains 'sp.'"
        if 'cf.' in sp_name.lower():
            return False, "name contains 'cf.'"
        if ' endosymbiont' in sp_name.lower():
            return False, "name contains 'endosymbiont'"
        if ' symbiont' in sp_name.lower():
            return False, "name contains 'symbiont'"
        if ' mycovirus' in sp_name.lower():
            return False, "name contains 'mycovirus'"
        if sp_name.lower().split()[1] == 'oral':
            return False, "specific name is 'oral'"
        if 'candidatus' in sp_name.lower() and sp_name.lower().split()[2] == 'oral':
            return False, "specific name is 'oral'"
        if '-like' in test_name.lower():
            return False, "full name contains '-like'"
        if 'endosymbiont' in test_name.lower().split():
            return False, "full name contains 'endosymbiont'"
        if 'symbiont' in test_name.lower().split():
            return False, "full name contains 'symbiont'"
        if 'mycovirus' in test_name.lower().split():
            return False, "full name contains 'mycovirus'"
        if 'phytoplasma' in test_name.split():
            # note the Phytoplasma is a valid genus so we are
            # specifically looking for a lowercase 'p'
            return False, "full name contains 'phytoplasma'"

        # check that binomial name contains only valid characters
        for ch in sp_name:  # ***
            if not ch.isalpha() and ch not in [' ', '[', ']']:
                return False, 'species name contains invalid character'

        return True, 's__' + sp_name

    def standardize_taxonomy(self, ncbi_taxonomy_file, output_consistent):
        """Produce standardized 7-rank taxonomy file from NCBI taxonomy strings."""

        fout_consistent = open(output_consistent, 'w')
        failed_filters = set()
        for line in open(ncbi_taxonomy_file):
            line_split = line.strip().split('\t')

            gid = line_split[0]
            taxonomy = line_split[1].split(';')

            if not ('d__Bacteria' in taxonomy or 'd__Archaea' in taxonomy):
                continue

            # remove unrecognized ranks (i.e., 'x__') and strain classification
            revised_taxonomy = []
            for t in taxonomy:
                if not t.startswith('x__') and not t.startswith('st__'):
                    revised_taxonomy.append(t)

            # create longest taxonomy string possible with canonical ranks
            canonical_taxonomy = {}
            for i, taxon in enumerate(revised_taxonomy):
                rank_prefix = taxon[0:3]
                if rank_prefix in Taxonomy.rank_prefixes:
                    if rank_prefix == 's__':
                        valid_name, canonical_species_name = self._valid_species_name(
                            taxon)

                        if valid_name:
                            canonical_taxonomy[Taxonomy.rank_prefixes.index(
                                rank_prefix)] = canonical_species_name
                        else:
                            if ('full name' in canonical_species_name and
                                ('oral' in canonical_species_name
                                 or '-like' in canonical_species_name
                                 or 'endosymbiont' in canonical_species_name
                                 or 'symbiont' in canonical_species_name
                                 or 'mycovirus' in canonical_species_name
                                 or 'phytoplasma' in canonical_species_name)):
                                failed_filters.add(taxon)
                    else:
                        canonical_taxonomy[Taxonomy.rank_prefixes.index(
                            rank_prefix)] = taxon

            # fill in missing ranks where possible
            if canonical_taxonomy:
                for i in xrange(0, max(canonical_taxonomy.keys())):
                    if i in canonical_taxonomy and (i + 1) not in canonical_taxonomy:
                        canonical_taxonomy[i +
                                           1] = Taxonomy.rank_prefixes[i + 1]

            cur_taxonomy = []
            for i in xrange(0, len(Taxonomy.rank_prefixes)):
                if i in canonical_taxonomy:
                    cur_taxonomy.append(canonical_taxonomy[i])
                else:
                    break  # unable to correctly determine a valid taxonomy below this rank

            if len(cur_taxonomy) > 0:
                if len(cur_taxonomy) != len(Taxonomy.rank_prefixes):
                    cur_taxonomy = cur_taxonomy + \
                        list(Taxonomy.rank_prefixes[len(cur_taxonomy):])
                fout_consistent.write('%s\t%s\n' %
                                      (gid, ';'.join(cur_taxonomy)))

        fout_consistent.close()

        # Sanity check particular filters
        fout = open('failed_filters.tsv', 'w')
        for sp in failed_filters:
            fout.write(sp + '\n')
        fout.close()

        print 'Genomes with a consistent taxonomy written to: %s' % output_consistent

    def run(self,
            taxonomy_dir,
            refseq_archaea_assembly_file,
            refseq_bacteria_assembly_file,
            genbank_archaea_assembly_file,
            genbank_bacteria_assembly_file,
            output_prefix):
        """Read NCBI taxonomy information and create summary output files."""

        # parse organism name
        self._assembly_organism_name(refseq_archaea_assembly_file,
                                     refseq_bacteria_assembly_file,
                                     genbank_archaea_assembly_file,
                                     genbank_bacteria_assembly_file,
                                     output_prefix + '_organism_names.tsv')

        # parse metadata file and taxonomy files
        assembly_to_tax_id = self._assembly_to_tax_id(refseq_archaea_assembly_file,
                                                      refseq_bacteria_assembly_file,
                                                      genbank_archaea_assembly_file,
                                                      genbank_bacteria_assembly_file)

        node_records = self._read_nodes(
            os.path.join(taxonomy_dir, 'nodes.dmp'))
        print 'Read %d node records.' % len(node_records)

        name_records = self._read_names(
            os.path.join(taxonomy_dir, 'names.dmp'))
        print 'Read %d name records.' % len(name_records)

        # traverse taxonomy tree for each assembly
        taxonomy_file = output_prefix + '_unfiltered_taxonomy.tsv'
        fout = open(taxonomy_file, 'w')

        print 'Number of assemblies: %d' % len(assembly_to_tax_id)
        for assembly_accession, tax_id in assembly_to_tax_id.iteritems():
            # traverse taxonomy tree to the root which is 'cellular organism' for genomes,
            # 'other sequences' for plasmids, and 'unclassified sequences' for metagenomic libraries
            taxonomy = []
            cur_tax_id = tax_id

            if cur_tax_id not in name_records:
                print '[Warning] Assembly %s has an invalid taxid: %s' % (assembly_accession, tax_id)
                continue

            roots = ['cellular organisms', 'other sequences',
                     'unclassified sequences', 'Viruses', 'Viroids']
            while name_records[cur_tax_id].name_txt not in roots:
                if cur_tax_id == '1':
                    print '[Error] TaxId %s reached root of taxonomy tree: %s' % (tax_id, taxonomy)
                    sys.exit(-1)

                try:
                    node_record = node_records[cur_tax_id]

                    if node_record.rank in Taxonomy.rank_labels:
                        rank_index = Taxonomy.rank_labels.index(
                            node_record.rank)
                        rank_prefix = Taxonomy.rank_prefixes[rank_index]
                    elif node_record.rank == 'subspecies':
                        rank_prefix = 'sb__'
                    else:
                        # unrecognized rank
                        rank_prefix = 'x__'
                        if node_record.rank == 'superkingdom':
                            rank_prefix = 'd__'

                    taxonomy.append(
                        rank_prefix + name_records[cur_tax_id].name_txt)

                    cur_tax_id = node_record.parent_tax_id
                except:
                    print traceback.format_exc()
                    print taxonomy

            taxonomy.reverse()
            taxa_str = ';'.join(taxonomy)
            fout.write('%s\t%s\n' % (assembly_accession, taxa_str))

        fout.close()

        self.standardize_taxonomy(taxonomy_file,
                                  output_prefix + '_standardized.tsv')


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'taxonomy_dir', help='directory containing NCBI taxonomy files')
    parser.add_argument('refseq_archaea_assembly_file',
                        help='file with metadata for each RefSeq assembly')
    parser.add_argument('refseq_bacteria_assembly_file',
                        help='file with metadata for each GenBank assembly')
    parser.add_argument('genbank_archaea_assembly_file',
                        help='file with metadata for each GenBank assembly')
    parser.add_argument('genbank_bacteria_assembly_file',
                        help='file with metadata for each GenBank assembly')
    parser.add_argument('output_prefix', help='output prefix')

    args = parser.parse_args()

    try:
        p = TaxonomyNCBI()
        p.run(args.taxonomy_dir,
              args.refseq_archaea_assembly_file,
              args.refseq_bacteria_assembly_file,
              args.genbank_archaea_assembly_file,
              args.genbank_bacteria_assembly_file,
              args.output_prefix)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
