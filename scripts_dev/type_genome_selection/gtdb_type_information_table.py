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

__prog_name__ = 'gtdb_type_information_table.py'
__prog_desc__ = 'Parse strain repositories to identify genomes assembled from type material.'

__author__ = 'Pierre Chaumeil and Donovan Parks'
__copyright__ = 'Copyright 2018'
__credits__ = ['Pierre Chaumeil', 'Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.5'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

""" VERSION HISTORY

0.0.5
- fixed issue with priority date not being set for type strain of species when
  the strain match was through an "unofficial" name

"""

import sys
import argparse
import os
import logging
import re
import time
from itertools import islice
from datetime import datetime
import math
import pickle
from collections import Counter, defaultdict, namedtuple
import multiprocessing as mp

from biolib.logger import logger_setup


class InfoGenerator(object):
    """Parse multiple sources to identify genomes assembled from type material."""

    def __init__(self, cpus, output_dir):
        """Initialization."""

        self.output_dir = output_dir
        self.cpus = cpus

        self.TYPE_SPECIES = 'type strain of species'
        self.TYPE_NEOTYPE = 'type strain of neotype'
        self.TYPE_SUBSPECIES = 'type strain of subspecies'
        self.TYPE_HETERO_SYNONYM = 'type strain of heterotypic synonym'
        self.NOT_TYPE_MATERIAL = 'not type material'
        self.type_priority = [self.TYPE_SPECIES,
                              self.TYPE_NEOTYPE,
                              self.TYPE_SUBSPECIES,
                              self.TYPE_HETERO_SYNONYM,
                              self.NOT_TYPE_MATERIAL]

        self.Match = namedtuple('Match', ['category',
                                          'istype',
                                          'isneotype',
                                          'gtdb_type_status',
                                          'standard_name',
                                          'strain_id',
                                          'year_date'])

        logger_setup(args.output_dir, "gtdb_type_information_table.log",
                     "gtdb_type_information_table", __version__, False)
        self.logger = logging.getLogger('timestamp')

    def load_year_dict(self, year_table):
        """Load year of priority for species as identified at LPSN, DSMZ, and StrainInfo."""

        dict_date = {}
        with open(year_table) as yt:
            for line in yt:
                infos = line.rstrip('\n').split('\t')
                dict_date[infos[0]] = {'lpsn': infos[1],
                                       'dsmz': infos[2],
                                       'straininfo': infos[3]}
        return dict_date

    def _standardize_strain_id(self, strain_id):
        """Convert strain ID into standard format."""

        pattern = re.compile('[\W_]+')
        strain_id = strain_id.replace('strain', '')
        standardized_id = pattern.sub('', strain_id.strip()).upper()

        return standardized_id

    def _fix_common_strain_id_errors(self, strain_ids):
        """Fix common erros associated with NCBI strain IDs."""

        # NCBI strain IDs sometimes contain a 'T' at the end that
        # actually designates the ID is for type material. To
        # resolve this tailing T's are removed and the new ID
        # added as a potential strain ID
        new_ids = set()
        for sid in strain_ids:
            if len(sid) > 1 and sid[-1] == 'T':
                new_ids.add(sid[0:-1])

            new_ids.add(sid)

        return new_ids

    def parse_ncbi_names_and_nodes(self, ncbi_names_file, ncbi_nodes_file, taxids_of_interest):
        """Parse NCBI names.dmp and nodes.dmp files"""

        # determine NCBI taxIDs of species and parent<->child tree
        species_taxids = set()
        parent = {}
        for line in open(ncbi_nodes_file):
            tokens = [token.strip() for token in line.split('|')]

            cur_taxid = int(tokens[0])
            parent_taxid = int(tokens[1])
            rank = tokens[2]

            parent[cur_taxid] = parent_taxid

            if rank == 'species':
                species_taxids.add(cur_taxid)

        self.logger.info(
            'Identified %d NCBI taxonomy species nodes.' % len(species_taxids))

        # determine species taxID of all taxa of interest
        species_of_taxid = {}
        for cur_taxid in taxids_of_interest:
            parent_taxid = cur_taxid
            while True:
                if parent_taxid in species_taxids:
                    species_of_taxid[cur_taxid] = parent_taxid
                    break

                if parent_taxid not in parent:
                    # this happens as not all genomes are defined below
                    # the rank of species and since the NCBI taxonomy and
                    # genome data are not always in sync
                    break

                parent_taxid = parent[parent_taxid]

        self.logger.info(
            'Associated %d NCBI taxon nodes with their parent species node.' % len(species_of_taxid))

        # parse auxillary names associated with a NCBI taxID and
        # type material strain IDs for species
        category_names = {}
        type_material = defaultdict(set)
        ncbi_authority = {}
        with open(ncbi_names_file) as nnf:
            for line in nnf:
                tokens = [token.strip() for token in line.split('|')]
                cur_taxid = int(tokens[0])

                if tokens[3] == 'type material':
                    for sid in self._fix_common_strain_id_errors([tokens[1]]):
                        type_material[cur_taxid].add(
                            self._standardize_strain_id(sid))

                if cur_taxid in taxids_of_interest:
                    if tokens[3] == 'authority':
                        ncbi_authority[cur_taxid] = tokens[1]
                    if tokens[3] in ['misspelling', 'synonym', 'equivalent name', 'scientific name']:
                        if cur_taxid not in category_names:
                            category_names[cur_taxid] = {'misspelling': [],
                                                         'synonym': [],
                                                         'equivalent name': [],
                                                         'scientific name': []}
                        category_names[cur_taxid][tokens[3]].append(tokens[1])

        self.logger.info(
            'Read auxillary species name information for %d NCBI taxIDs.' % len(category_names))
        self.logger.info(
            'Read type material information for %d NCBI taxIDs.' % len(type_material))

        # sanity check results
        for k, v in category_names.iteritems():
            if len(set(v['synonym']).intersection(v.get('scientific name'))) > 0 or len(set(v['synonym']).intersection(v['equivalent name'])) > 0:
                print 'ERROR'
                print v['synonym']
                print v.get('scientific name')
                print v['equivalent name']
                sys.exit(-1)

        return category_names, type_material, species_of_taxid, ncbi_authority

    def load_metadata(self, metadata_file):
        """Parse data from GTDB metadata file."""

        metadata = {}
        with open(metadata_file) as metaf:
            headers_line = metaf.readline()
            separator = ','
            if '\t' in headers_line:
                separator = '\t'

            headers = headers_line.rstrip('\n').split(separator)

            gtdb_ncbi_organism_name_index = headers.index('ncbi_organism_name')
            gtdb_ncbi_type_material_designation_index = headers.index(
                'ncbi_type_material_designation')
            gtdb_accession_index = headers.index('accession')
            gtdb_taxonomy_species_name_index = headers.index('ncbi_taxonomy')
            gtdb_strain_identifiers_index = headers.index(
                'ncbi_strain_identifiers')
            gtdb_ncbi_taxonomy_unfiltered_index = headers.index(
                'ncbi_taxonomy_unfiltered')
            gtdb_ncbi_taxid_index = headers.index('ncbi_taxid')

            taxids = set()
            for line in metaf:
                infos = line.rstrip('\n').split(separator)

                if not infos[gtdb_accession_index].startswith('U_'):
                    # standardize NCBI strain IDs
                    standard_strain_ids = []
                    if infos[gtdb_strain_identifiers_index] != 'none':
                        pattern = re.compile('[\W_]+')
                        created_list = [
                            sid.strip() for sid in infos[gtdb_strain_identifiers_index].split(';')]
                        created_list = self._fix_common_strain_id_errors(
                            created_list)
                        standard_strain_ids = [self._standardize_strain_id(sid)
                                               for sid in created_list
                                               if (sid != '' and sid != 'none')]

                    metadata[infos[gtdb_accession_index]] = {
                        'ncbi_organism_name': infos[gtdb_ncbi_organism_name_index],
                        'taxonomy_species_name': infos[gtdb_taxonomy_species_name_index].split(';')[6].replace('s__',
                                                                                                               ''),
                        'ncbi_strain_ids': infos[gtdb_strain_identifiers_index],
                        'ncbi_standardized_strain_ids': set(standard_strain_ids),
                        'ncbi_type_material_designation': infos[gtdb_ncbi_type_material_designation_index],
                        'ncbi_taxonomy_unfiltered': infos[gtdb_ncbi_taxonomy_unfiltered_index],
                        'ncbi_taxid': int(infos[gtdb_ncbi_taxid_index])}

                    taxids.add(int(infos[gtdb_ncbi_taxid_index]))

        return metadata, taxids

    def load_straininfo_strains_dictionary(self, straininfo_dir):
        # We load the dictionary of strains from Straininfo
        straininfo_strains_dic = {}
        pattern = re.compile('[\W_]+')
        p = re.compile('(\w+)\s\(now\s(\w+)\)\s(\d+)', re.IGNORECASE)
        with open(os.path.join(straininfo_dir, 'straininfo_strains.tsv')) as ststr:
            ststr.readline()
            for line in ststr:
                infos = line.rstrip('\n').split('\t')
                if len(infos) < 2:
                    print "len(infos) < 2 "
                    print infos
                else:
                    new_ls = []
                    list_strains = infos[1].split("=")
                    for ite in list_strains:
                        new_ls.append(ite)
                        matches = p.search(ite)
                        # if the strain is similar to "IFO (now NBC) 12345"
                        # we store the strains IFO12345, IFO 12345, NBC 12345, NCB12345
                        # straininfo is the only souces where we have this format
                        # dsmz and lpsn are already processed
                        if matches:
                            new_ls.append("{} {}".format(
                                matches.group(1), matches.group(3)))
                            new_ls.append("{}{}".format(
                                matches.group(1), matches.group(3)))
                            new_ls.append("{} {}".format(
                                matches.group(2), matches.group(3)))
                            new_ls.append("{}{}".format(
                                matches.group(2), matches.group(3)))

                        standard_list_strains = [pattern.sub('', a.strip()).upper(
                        ) for a in new_ls if (a != '' and a != 'none')]

                        straininfo_strains_dic[infos[0]] = "=".join(
                            set(standard_list_strains))
        return straininfo_strains_dic

    def load_dsmz_strains_dictionary(self, dsmz_dir):
        # We load the dictionary of strains from DSMZ
        dsmz_strains_dic = {}
        pattern = re.compile('[\W_]+')
        with open(os.path.join(dsmz_dir, 'dsmz_strains.tsv')) as dsstr:
            dsstr.readline()
            for line in dsstr:
                infos = line.rstrip('\n').split('\t')
                if len(infos) < 2:
                    print "len(infos) < 2 "
                    print infos
                else:
                    list_strains = [pattern.sub('', a.strip()).upper(
                    ) for a in infos[1].split('=') if (a != '' and a != 'none')]
                    dsmz_strains_dic[infos[0]] = '='.join(set(list_strains))

        return dsmz_strains_dic

    def load_lpsn_strains_dictionary(self, lpsn_dir):
        # We load the dictionary of strains from LPSN
        pattern = re.compile('[\W_]+')
        lpsn_strains_dic = {}
        with open(os.path.join(lpsn_dir, 'lpsn_strains.tsv')) as lpstr:
            lpstr.readline()
            for line in lpstr:
                infos = line.rstrip('\n').split('\t')

                if len(infos) < 3:
                    print "len(infos) < 3 "
                    print infos
                else:
                    list_strains = [pattern.sub('', a.strip()).upper(
                    ) for a in infos[1].split('=') if (a != '' and a != 'none')]
                    list_neotypes = [pattern.sub('', a.strip()).upper(
                    ) for a in infos[2].split('=') if (a != '' and a != 'none')]

                    lpsn_strains_dic[infos[0]] = {
                        'strains': '='.join(set(list_strains)), 'neotypes': '='.join(set(list_neotypes))}
        return lpsn_strains_dic

    def _read_type_species_of_genus(self, species_file):
        """Read type species of genus."""

        type_species_of_genus = {}
        genus_type_species = {}
        with open(species_file) as lpstr:
            lpstr.readline()

            for line in lpstr:
                line_split = line.rstrip('\n').split('\t')
                sp, genus, authority = line_split

                if genus:
                    sp = sp.replace('s__', '')
                    type_species_of_genus[sp] = genus

                    if genus in genus_type_species:
                        self.logger.warning(
                            'Identified multiple type species for %s in %s.' % (genus, species_file))
                    else:
                        genus_type_species[genus] = sp

        return type_species_of_genus, genus_type_species

    def remove_brackets(self, sp_name):
        """Remove brackets from species name.

        e.g., st__[Eubacterium] siraeum 70/3
        """

        if sp_name.startswith('['):
            sp_name = sp_name.replace('[', '', 1).replace(']', '', 1)
        return sp_name

    def get_species_name(self, gid):
        """Determine species name for genome.

        This is the NCBI subspecies name if defined,
        and the species name otherwise.
        """

        ncbi_unfiltered_taxa = [
            t.strip() for t in self.metadata[gid]['ncbi_taxonomy_unfiltered'].split(';')]
        ncbi_species = None
        ncbi_subspecies = None
        for taxon in ncbi_unfiltered_taxa:
            if taxon.startswith('s__'):
                ncbi_species = taxon[3:]
            elif taxon.startswith('sb__'):
                ncbi_subspecies = taxon[4:]

                # fix odd designation impacting less than a dozen genomes
                ncbi_subspecies = ncbi_subspecies.replace(' pv. ', ' subsp. ')

        if ncbi_subspecies:
            if 'subsp.' not in ncbi_subspecies:
                self.logger.warning(
                    "NCBI subspecies name without 'subsp.' definition: %s" % ncbi_subspecies)

            return self.remove_brackets(ncbi_subspecies)

        if ncbi_species:
            return self.remove_brackets(ncbi_species)

        return None

    def get_priority_year(self, spe_name, sourcest):
        """Get year of priority for species."""

        if spe_name in self.year_table:
            if self.year_table[spe_name][sourcest] != '':
                return int(self.year_table[spe_name][sourcest])

        return ''

    def strains_iterate(self, gid, standard_name, repository_strain_ids, raw_names, misspelling_names, synonyms, equivalent_names, isofficial, sourcest):
        """Check for matching species name and type strain IDs."""

        # search each strain ID at a given strain repository (e.g. LPSN)
        # associated with the species name
        istype = False
        year_date = ''
        category_name = ''
        matched_strain_id = None
        for repository_strain_id in repository_strain_ids.split("="):
            strain_ids = self.metadata[gid]['ncbi_expanded_standardized_strain_ids']
            if repository_strain_id in strain_ids:
                istype = True
            else:
                if len(repository_strain_id) <= 1:
                    continue  # too short to robustly identify

                # remove all white spaces and underscores, and capitalize, before
                # looking for match with standardized strain ID
                pattern = re.compile('[\W_]+')
                collapsed_names = {pattern.sub(
                    '', a).upper(): a for a in raw_names}

                for name in collapsed_names:
                    if repository_strain_id in name:
                        first_char = repository_strain_id[0]
                        p = re.compile(' {}'.format(first_char), re.IGNORECASE)
                        matches_beginning = p.search(collapsed_names[name])

                        last_char = repository_strain_id[-1]
                        q = re.compile('{}(\s|$)'.format(
                            last_char), re.IGNORECASE)
                        matches_end = q.search(collapsed_names[name])
                        if matches_beginning and matches_end:
                            istype = True

            if istype:
                year_date = self.get_priority_year(standard_name, sourcest)
                
                if not isofficial:
                    category_name = self.select_category_name(standard_name,
                                                              misspelling_names,
                                                              synonyms,
                                                              equivalent_names)
                else:
                    category_name = 'official_name'

                matched_strain_id = repository_strain_id
                break

        return (matched_strain_id, category_name, istype, year_date)

    def type_species_or_subspecies(self, gid):
        """Determine if genome is the 'type strain of species' or 'type strain of subspecies'."""

        sp_name = self.get_species_name(gid)
        if 'subsp.' not in sp_name:
            return self.TYPE_SPECIES
        else:
            tokens = sp_name.split()
            subsp_index = tokens.index('subsp.')
            if tokens[subsp_index - 1] == tokens[subsp_index + 1]:
                return self.TYPE_SPECIES

        return self.TYPE_SUBSPECIES

    def match_with_latinization(self, test_sp, target_sp):
        """Check for a match between the specific name of two species considering different gender suffixes."""

        # get specific name from species name
        test = test_sp.split()[1]
        target = target_sp.split()[1]
        if test == target:
            return True

        # determine gender of test name and check for match
        # with related suffix from same group of Latin adjectives
        masc = ('us', 'is', 'er')
        fem = ('a', 'is', 'eris')
        neu = ('um', 'e', 'ere')
        for s, s1, s2 in [(masc, fem, neu), (fem, masc, neu), (neu, masc, fem)]:
            for idx, suffix in enumerate(s):
                if test.endswith(suffix):
                    if test[0:-len(suffix)] + s1[idx] == target:
                        return True
                    elif test[0:-len(suffix)] + s2[idx] == target:
                        return True

        return False

    def check_heterotypic_synonym(self, spe_name, official_spe_names):
        """Check if species is a heterotypic synonym (i.e. has difference specific name)."""

        for official_name in official_spe_names:
            if self.match_with_latinization(spe_name, official_name):
                return False

        return True

    def strain_match(self,
                     gid,
                     standard_names,
                     official_standard_names,
                     misspelling_names,
                     synonyms,
                     equivalent_names,
                     strain_dictionary,
                     sourcest,
                     isofficial):
        """Match species names with stain IDs for a type source (e.g. LPSN) in order to establish if a genome is assembled from type."""

        # Match strain IDs from type sources (e.g. LPSN) associated with each standard
        # species name to strain information at NCBI. Searching is performed on the
        # raw NCBI species designations associated with a standard name as it can be
        # challenging to parse strain information from these entries.
        match = None
        gtdb_types = set()
        for standard_name, raw_names in standard_names.items():
            if standard_name not in strain_dictionary:
                continue

            if sourcest == 'lpsn':
                # lpsn has information for both strains and neotype strains
                repository_strain_ids = strain_dictionary.get(
                    standard_name).get('strains')
            else:
                repository_strain_ids = strain_dictionary.get(standard_name)

            matched_strain_id, category, istype, year_date = self.strains_iterate(gid, 
                                                                standard_name, 
                                                                repository_strain_ids, 
                                                                raw_names, 
                                                                misspelling_names, 
                                                                synonyms, 
                                                                equivalent_names, 
                                                                isofficial, 
                                                                sourcest)
            

            isneotype = False
            if not istype and sourcest == 'lpsn':
                repository_strain_ids = strain_dictionary.get(
                    standard_name).get('neotypes')
                matched_strain_id, _, isneotype, _ = self.strains_iterate(gid,
                                                                          standard_name,
                                                                          repository_strain_ids,
                                                                          raw_names,
                                                                          misspelling_names,
                                                                          synonyms,
                                                                          equivalent_names,
                                                                          isofficial,
                                                                          sourcest)

            gtdb_type_status = 'not type material'
            if istype or isneotype:
                heterotypic_synonym = False
                if not isofficial:
                    heterotypic_synonym = self.check_heterotypic_synonym(
                        standard_name, official_standard_names)

                if heterotypic_synonym:
                    gtdb_type_status = self.TYPE_HETERO_SYNONYM
                else:
                    gtdb_type_status = self.type_species_or_subspecies(gid)
                    if isneotype and gtdb_type_status == self.TYPE_SPECIES:
                        gtdb_type_status = self.TYPE_NEOTYPE

            m = self.Match(category, istype, isneotype, gtdb_type_status,
                           standard_name, matched_strain_id, year_date)
            if category == 'official_name':
                if match:
                    prev_gtdb_type_status = match.gtdb_type_status
                    if prev_gtdb_type_status != gtdb_type_status:
                        self.logger.error('Official species name has ambiguous type status: %s, %s, %s' % (
                            gid,
                            gtdb_type_status,
                            prev_gtdb_type_status))
                        sys.exit(-1)
                match = m
            elif category != '':
                # it is possible for a genome to be both a 'type strain of subspecies',
                # 'type strain of heterotypic synonym', and potentially a 'type strain of species'
                # depending on the different synonyms, equivalent names, and
                # misspelling
                if match:
                    prev_gtdb_type_status = match.gtdb_type_status
                    if self.type_priority.index(gtdb_type_status) < self.type_priority.index(prev_gtdb_type_status):
                        match = m
                else:
                    match = m

        return match

    def select_category_name(self, spe_name, misspelling_names, synonyms, equivalent_names):
        """Determine if name is a synonym, equivalent name, or misspelling."""

        # determine source of name giving highest priority to synonyms
        # and lowest priority to misspelling
        if spe_name in synonyms:
            return 'synonyms'
        elif spe_name in equivalent_names:
            return 'equivalent name'
        elif spe_name in misspelling_names:
            return 'misspelling name'

        self.logger.error('Failed to identify category of name: %s' % spe_name)
        sys.exit(-1)

    def standardise_names(self, potential_names):
        """Create a standard set of species names, include subsp. designations."""

        standardized = defaultdict(set)
        for raw_name in potential_names:

            # standardize the species name
            standard_name = re.sub(r'(?i)(candidatus\s)', r'', raw_name)
            standard_name = re.sub(r'(?i)(serotype.*)', r'', standard_name)
            standard_name = re.sub(r'(?i)(ser\..*)', r'', standard_name)
            standard_name = re.sub(r'\"|\'', r'', standard_name)
            standard_name = self.remove_brackets(standard_name)
            standard_name = standard_name.strip()

            name_tokens = standard_name.split(' ')

            # check if name is binomial
            if len(name_tokens) != 2:
                if len(name_tokens) == 4 and name_tokens[2] == 'subsp.':
                    standardized[standard_name].add(raw_name)

                    # if the subspecies matches the species name
                    if name_tokens[1] == name_tokens[3]:
                        standardized[' '.join(name_tokens[0:2])].add(raw_name)

                # if the name is longer than 4 words but the 3rd word still
                # subsp, we assume than the 4 first words are a subspecies name
                elif len(name_tokens) >= 4 and name_tokens[2] == 'subsp.':
                    if name_tokens[1] == name_tokens[3]:
                        standardized[' '.join(name_tokens[0:2])].add(raw_name)
                        subsp_name = '{0} {1} subsp. {1}'.format(name_tokens[0],
                                                                 name_tokens[1])
                        standardized[subsp_name].add(raw_name)
                    if name_tokens[1] != name_tokens[3]:
                        standardized[' '.join(name_tokens[0:4])].add(raw_name)
                elif len(name_tokens) >= 2:
                    standardized[' '.join(name_tokens[0:2])].add(raw_name)
                    subsp_name = '{0} {1} subsp. {1}'.format(name_tokens[0],
                                                             name_tokens[1])
                    standardized[subsp_name].add(raw_name)
            else:
                standardized[standard_name].add(standard_name)
                subsp_name = '{0} {1} subsp. {1}'.format(name_tokens[0],
                                                         name_tokens[1])
                standardized[subsp_name].add(raw_name)

        return standardized

    def parse_strains(self, sourcest, strain_dictionary, outfile):
<<<<<<< Updated upstream
        """Parse information for a single strain resource (e.g., LPSN, DSMZ, or StrainInfo)."""
        
=======
        """Parse information for a single strain resouce (e.g., LPSN, DSMZ, or StrainInfo)."""

>>>>>>> Stashed changes
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for gid in self.metadata:
            worker_queue.put(gid)

        for _ in range(self.cpus):
            worker_queue.put(None)

        try:
            workerProc = [mp.Process(target=self._worker, args=(sourcest,
                                                                strain_dictionary,
                                                                worker_queue,
                                                                writer_queue)) for _ in range(self.cpus)]
            writeProc = mp.Process(target=self._writer, args=(
                sourcest, outfile, writer_queue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writer_queue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()
            writeProc.terminate()

    def _worker(self,
                sourcest,
                strain_dictionary,
                queue_in,
                queue_out):
        """Determine if genome is assembled from type material."""

        while True:
            gid = queue_in.get(block=True, timeout=None)
            if gid == None:
                break

            genome_metadata = self.metadata[gid]

            species_name = self.get_species_name(gid)
            if species_name is None:
                continue
            standardized_sp_names = self.standardise_names([species_name])

            # get list of misspellings, synonyms, and equivalent names associated
            # with this genome
            unofficial_potential_names = set()
            if genome_metadata['ncbi_taxid'] in self.ncbi_auxiliary_names:
                unofficial_potential_names.update(self.ncbi_auxiliary_names[
                    genome_metadata['ncbi_taxid']]['misspelling'])
                unofficial_potential_names.update(self.ncbi_auxiliary_names[
                    genome_metadata['ncbi_taxid']]['synonym'])
                unofficial_potential_names.update(self.ncbi_auxiliary_names[
                    genome_metadata['ncbi_taxid']]['equivalent name'])

                misspelling_names = self.standardise_names(self.ncbi_auxiliary_names[
                    genome_metadata['ncbi_taxid']]['misspelling'])
                synonyms = self.standardise_names(self.ncbi_auxiliary_names[
                    genome_metadata['ncbi_taxid']]['synonym'])
                equivalent_names = self.standardise_names(self.ncbi_auxiliary_names[
                    genome_metadata['ncbi_taxid']]['equivalent name'])

            unofficial_standard_names = self.standardise_names(
                unofficial_potential_names)

            # match species and strain information from NCBI with information
            # at type repository (e.g., LPSN)
            match = self.strain_match(gid,
                                      standardized_sp_names,
                                      standardized_sp_names,
                                      None,
                                      None,
                                      None,
                                      strain_dictionary,
                                      sourcest,
                                      True)

            if not match:
                # check if any of the auxillary names have a species name
                # and strain ID match with the type repository
                match = self.strain_match(gid,
                                          unofficial_standard_names,
                                          standardized_sp_names,
                                          misspelling_names,
                                          synonyms,
                                          equivalent_names,
                                          strain_dictionary,
                                          sourcest,
                                          False)

            if match:
                if sourcest == 'lpsn':
                    # lpsn has information for both strains and neotype strains
                    repository_strain_ids = strain_dictionary[match.standard_name].get(
                        'strains')
                else:
                    repository_strain_ids = strain_dictionary[match.standard_name]

                queue_out.put((gid,
                               species_name,
                               match.year_date,
                               match.istype,
                               match.isneotype,
                               match.gtdb_type_status,
                               match.category,
                               match.standard_name,
                               match.strain_id,
                               set(repository_strain_ids.split('='))))

    def _writer(self, sourcest, outfile, writer_queue):
        """Report type material status for each genome."""

        fout = open(outfile, 'w')
        fout.write(
            'genome\tncbi_organism_name\tncbi_species_name\tncbi_type_designation\tgtdb_type_designation')
        fout.write(
            '\tncbi_base_strain_ids\tncbi_canonical_strain_ids\tmatched_strain_id')
        fout.write(
            '\t{0}_match_type\t{0}_match_name\t{0}_match_strain_id\t{0}_strain_ids'.format(sourcest))
        fout.write('\tmissspellings\tequivalent_names\tsynonyms')
        fout.write('\tneotype\tpriority_year\n')

        processed = 0
        while True:
            data = writer_queue.get(block=True, timeout=None)
            if data == None:
                break

            (gid,
                species_name,
                year_date,
                type_strain,
                neotype,
                gtdb_type_status,
                category_name,
                matched_sp_name,
                matched_strain_id,
                repository_strain_ids,) = data

            info_genomes = self.metadata[gid]

            misspelling = equivalent_name = synonym = ''
            if info_genomes['ncbi_taxid'] in self.ncbi_auxiliary_names:
                misspelling = '; '.join(
                    self.ncbi_auxiliary_names[info_genomes['ncbi_taxid']]['misspelling'])
                equivalent_name = '; '.join(
                    self.ncbi_auxiliary_names[info_genomes['ncbi_taxid']]['equivalent name'])
                synonym = '; '.join(
                    self.ncbi_auxiliary_names[info_genomes['ncbi_taxid']]['synonym'])

            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                gid,
                info_genomes['ncbi_organism_name'],
                species_name,
                info_genomes['ncbi_type_material_designation'],
                gtdb_type_status,
                self.metadata[gid]['ncbi_strain_ids'],
                '; '.join(
                    self.metadata[gid]['ncbi_expanded_standardized_strain_ids']),
                matched_strain_id,
                category_name,
                matched_sp_name,
                '; '.join(repository_strain_ids.intersection(
                    self.metadata[gid]['ncbi_expanded_standardized_strain_ids'])),
                '; '.join(repository_strain_ids),
                misspelling,
                equivalent_name,
                synonym,
                neotype,
                year_date))

            processed += 1
            statusStr = '-> Processing %d genomes assembled from type material.'.ljust(
                86) % processed
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def _parse_strain_summary(self, strain_summary_file):
        """Parse type information from strain repository."""

        StrainInfo = namedtuple('StrainInfo', 'type_designation priority_year')

        strain_info = {}
        with open(strain_summary_file) as f:
            header = f.readline().rstrip().split('\t')

            gid_index = header.index('genome')
            gtdb_type_designation_index = header.index('gtdb_type_designation')
            priority_year_index = header.index('priority_year')

            for line in f:
                line_split = line.rstrip('\n').split('\t')

                gid = line_split[gid_index]
                type_designation = line_split[gtdb_type_designation_index]
                priority_year = line_split[priority_year_index]

                strain_info[gid] = StrainInfo(type_designation, priority_year)

        return strain_info

    def type_summary_table(self,
                           ncbi_authority,
                           lpsn_summary_file,
                           dsmz_summary_file,
                           straininfo_summary_file,
                           lpsn_type_species_of_genus,
                           dsmz_type_species_of_genus,
                           summary_table_file):
        """Generate type strain summary file across all strain repositories."""

        # parse strain repository files
        lpsn = self._parse_strain_summary(lpsn_summary_file)
        dsmz = self._parse_strain_summary(dsmz_summary_file)
        straininfo = self._parse_strain_summary(straininfo_summary_file)

        # write out type strain information for each genome
        fout = open(summary_table_file, 'w')
        fout.write(
            "accession\tncbi_species\tncbi_organism_name\tncbi_strain_ids\tncbi_canonical_strain_ids")
        fout.write("\tncbi_taxon_authority\tncbi_type_designation")
        fout.write("\tgtdb_type_designation\tgtdb_type_designation_sources")
        fout.write(
            "\tlpsn_type_designation\tdsmz_type_designation\tstraininfo_type_designation")
        fout.write(
            '\tlpsn_priority_year\tdsmz_priority_year\tstraininfo_priority_year')
        fout.write("\tgtdb_type_species_of_genus\n")

        missing_type_at_ncbi = 0
        missing_type_at_gtdb = 0
        agreed_type_of_species = 0
        agreed_type_of_subspecies = 0
        num_type_species_of_genus = 0
        for gid, metadata in self.metadata.iteritems():
            fout.write(gid)

            species_name = self.get_species_name(gid)
            fout.write('\t%s\t%s\t%s\t%s' % (species_name,
                                             metadata['ncbi_organism_name'],
                                             metadata['ncbi_strain_ids'],
                                             '; '.join(metadata['ncbi_expanded_standardized_strain_ids'])))

            fout.write('\t%s\t%s' % (ncbi_authority.get(metadata['ncbi_taxid'], '').replace('"', '~'),
                                     metadata['ncbi_type_material_designation']))

            # GTDB sets the type material designation in a specific priority
            # order
            highest_priority_designation = self.NOT_TYPE_MATERIAL
            for sr in [lpsn, dsmz, straininfo]:
                if gid in sr and self.type_priority.index(sr[gid].type_designation) < self.type_priority.index(highest_priority_designation):
                    highest_priority_designation = sr[gid].type_designation
            fout.write('\t%s' % highest_priority_designation)

            type_species_of_genus = False
            canonical_sp_name = ' '.join(species_name.split()[0:2])
            if (highest_priority_designation == 'type strain of species' and
                    (species_name in lpsn_type_species_of_genus or species_name in dsmz_type_species_of_genus
                        or canonical_sp_name in lpsn_type_species_of_genus or canonical_sp_name in dsmz_type_species_of_genus)):
                type_species_of_genus = True
                num_type_species_of_genus += 1

            gtdb_type_sources = []
            for sr_id, sr in [('LPSN', lpsn), ('DSMZ', dsmz), ('StrainInfo', straininfo)]:
                if gid in sr and sr[gid].type_designation == highest_priority_designation:
                    gtdb_type_sources.append(sr_id)
            fout.write('\t%s' % '; '.join(gtdb_type_sources))

            fout.write('\t%s\t%s\t%s' % (lpsn[gid].type_designation if gid in lpsn else self.NOT_TYPE_MATERIAL,
                                         dsmz[gid].type_designation if gid in dsmz else self.NOT_TYPE_MATERIAL,
                                         straininfo[gid].type_designation if gid in straininfo else self.NOT_TYPE_MATERIAL))
            fout.write('\t%s\t%s\t%s' % (lpsn[gid].priority_year if gid in lpsn else '',
                                         dsmz[gid].priority_year if gid in dsmz else '',
                                         straininfo[gid].priority_year if gid in straininfo else ''))
            fout.write('\t%s\n' % type_species_of_genus)

            if metadata['ncbi_type_material_designation'] == 'none' and highest_priority_designation == self.TYPE_SPECIES:
                missing_type_at_ncbi += 1

            if metadata['ncbi_type_material_designation'] == 'assembly from type material' and highest_priority_designation == self.NOT_TYPE_MATERIAL:
                missing_type_at_gtdb += 1

            if metadata['ncbi_type_material_designation'] == 'assembly from type material' and highest_priority_designation == self.TYPE_SPECIES:
                agreed_type_of_species += 1

            if (metadata['ncbi_type_material_designation'] in ['assembly from type material', 'assembly from synonym type material']
                    and highest_priority_designation == self.TYPE_SUBSPECIES):
                sp_tokens = self.get_species_name(gid).split()
                if sp_tokens[2] == 'subsp.' and sp_tokens[1] != sp_tokens[3]:
                    agreed_type_of_subspecies += 1

        self.logger.info(
            'Identified %d genomes designated as the type species of genus.' % num_type_species_of_genus)
        self.logger.info(
            'Genomes that appear to have missing type species information at NCBI: %d' % missing_type_at_ncbi)
        self.logger.info(
            'Genomes that are only effectively published or erroneously missing type species information at GTDB: %d' % missing_type_at_gtdb)
        self.logger.info(
            'Genomes where GTDB and NCBI both designate type strain of species: %d' % agreed_type_of_species)
        self.logger.info(
            'Genomes where GTDB and NCBI both designate type strain of subspecies: %d' % agreed_type_of_subspecies)

        fout.close()

    def _expand_ncbi_strain_ids(self, ncbi_coidentical_strain_ids, ncbi_species_of_taxid):
        """Expand set of NCBI co-identical strain IDs associated with each genome."""

        for gid, genome_metadata in self.metadata.items():
            # determine the list of strain IDs at NCBI that are
            # associated with the genome
            strain_ids = genome_metadata['ncbi_standardized_strain_ids']
            ncbi_taxid = genome_metadata['ncbi_taxid']
            if ncbi_taxid in ncbi_coidentical_strain_ids:
                if strain_ids.intersection(ncbi_coidentical_strain_ids[ncbi_taxid]):
                    # expand list of strain IDs to include all co-identical
                    # type material strain IDs specified by the NCBI taxonomy
                    # in names.dmp for this taxon
                    strain_ids = strain_ids.union(
                        ncbi_coidentical_strain_ids[ncbi_taxid])

            # check if genome is associated with a NCBI species node which may have
            # additional relevant co-identical strain IDs
            if ncbi_taxid in ncbi_species_of_taxid:
                ncbi_sp_taxid = ncbi_species_of_taxid[ncbi_taxid]
                if ncbi_sp_taxid in ncbi_coidentical_strain_ids:
                    if strain_ids.intersection(ncbi_coidentical_strain_ids[ncbi_sp_taxid]):
                        # expand list of strain IDs to include all co-identical
                        # type material strain IDs specified by the NCBI taxonomy
                        # in names.dmp for this taxon
                        strain_ids = strain_ids.union(
                            ncbi_coidentical_strain_ids[ncbi_sp_taxid])

            self.metadata[gid]['ncbi_expanded_standardized_strain_ids'] = strain_ids

    def run(self,
            metadata_file,
            ncbi_names_file,
            ncbi_nodes_file,
            lpsn_dir,
            dsmz_dir,
            straininfo_dir,
            year_table,
            sourcest):
        """Parse multiple sources to identify genomes assembled from type material."""

        # initialize data being parsed from file
        self.logger.info('Parsing GTDB metadata.')
        self.metadata, taxids_of_interest = self.load_metadata(metadata_file)

        self.logger.info('Parsing year table.')
        self.year_table = self.load_year_dict(year_table)

        self.logger.info(
            'Parsing NCBI taxonomy information from names.dmp and nodes.dmp.')
        rtn = self.parse_ncbi_names_and_nodes(
            ncbi_names_file, ncbi_nodes_file, taxids_of_interest)
        (self.ncbi_auxiliary_names,
            ncbi_coidentical_strain_ids,
            ncbi_species_of_taxid,
            ncbi_authority) = rtn

        # expand set of NCBI co-identical strain IDs associated with each
        # genome
        self.logger.info(
            'Expanding co-indetical strain IDs associated with each genome.')
        self._expand_ncbi_strain_ids(
            ncbi_coidentical_strain_ids, ncbi_species_of_taxid)

        # identify genomes assembled from type material
        self.logger.info('Identifying genomes assembled from type material.')
        if sourcest == 'lpsn' or sourcest == 'all':
            self.logger.info('Parsing information in LPSN directory.')
            self.lpsn_strains_dic = self.load_lpsn_strains_dictionary(lpsn_dir)

            self.logger.info('Processing LPSN data.')
            lpsn_summary_file = os.path.join(
                self.output_dir, 'lpsn_summary.tsv')
            self.parse_strains('lpsn',
                               self.lpsn_strains_dic,
                               lpsn_summary_file)

        if sourcest == 'dsmz' or sourcest == 'all':
            self.logger.info('Parsing information in DSMZ directory.')
            self.dsmz_strains_dic = self.load_dsmz_strains_dictionary(dsmz_dir)

            self.logger.info('Processing DSMZ data.')
            dsmz_summary_file = os.path.join(
                self.output_dir, 'dsmz_summary.tsv')
            self.parse_strains('dsmz',
                               self.dsmz_strains_dic,
                               dsmz_summary_file)

        if sourcest == 'straininfo' or sourcest == 'all':
            self.logger.info('Parsing information in StrainInfo directory.')
            self.straininfo_strains_dic = self.load_straininfo_strains_dictionary(
                straininfo_dir)

            self.logger.info('Processing StrainInfo data.')
            straininfo_summary_file = os.path.join(
                self.output_dir, 'straininfo_summary.tsv')
            self.parse_strains('straininfo',
                               self.straininfo_strains_dic,
                               straininfo_summary_file)

        # generate global summary file if information was generated from all
        # sources
        if sourcest == 'all':
            lpsn_type_species_of_genus, lpsn_genus_type_species = self._read_type_species_of_genus(
                os.path.join(lpsn_dir, 'lpsn_species.tsv'))
            dsmz_type_species_of_genus, dsmz_genus_type_species = self._read_type_species_of_genus(
                os.path.join(dsmz_dir, 'dsmz_species.tsv'))

            for genus in lpsn_genus_type_species:
                if genus in dsmz_genus_type_species:
                    if lpsn_genus_type_species[genus] != dsmz_genus_type_species[genus]:
                        self.logger.warning('LPSN and DSMZ disagree of type species for %s: %s %s' % (
                            genus,
                            lpsn_genus_type_species[genus],
                            dsmz_genus_type_species[genus]))

            self.logger.info('Identified %d LPSN and %d DSMZ type species of genus.' % (
                len(lpsn_type_species_of_genus),
                len(dsmz_type_species_of_genus)))

            self.logger.info(
                'Generating summary type information table across all strain repositories.')
            summary_table_file = os.path.join(
                self.output_dir, 'gtdb_type_strain_summary.tsv')
            self.type_summary_table(ncbi_authority,
                                    lpsn_summary_file,
                                    dsmz_summary_file,
                                    straininfo_summary_file,
                                    lpsn_type_species_of_genus,
                                    dsmz_type_species_of_genus,
                                    summary_table_file)

        self.logger.info('Done.')


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__
    print '  Contact: ' + __email__ + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--lpsn_dir',
                        help='Directory including the 3 LPSN result files (lpsn_genera.tsv, lpsn_species.tsv and lpsn_strains.tsv ).')
    parser.add_argument('--dsmz_dir',
                        help='Directory including the 3 DSMZ result files (dsmz_strains.tsv, dsmz_species.tsv and dsmz_genera.tsv ).')
    parser.add_argument('--straininfo_dir',
                        help='Directory including the straininfo result file (straininfo_strains.tsv).')
    parser.add_argument('--year_table',
                        help='Date table generated by generate_date_table.py.')
    parser.add_argument('--metadata_file',
                        help='Metadata file generated by GTDB')
    parser.add_argument('--ncbi_names', help='NCBI names.dmp file')
    parser.add_argument('--ncbi_nodes', help='NCBI nodes.dmp file')
    parser.add_argument('--cpus', help='Number of threads.',
                        type=int, default=1)
    parser.add_argument('--output_dir', help='Output directory.', default='.')
    parser.add_argument('--source_strain', choices=['all', 'lpsn', 'dsmz', 'straininfo'], default='all',
                        help='select LPSN,Straininfo,DSMZ to parse')

    args = parser.parse_args()

    try:
        p = InfoGenerator(args.cpus,
                          args.output_dir)

        p.run(args.metadata_file,
              args.ncbi_names,
              args.ncbi_nodes,
              args.lpsn_dir,
              args.dsmz_dir,
              args.straininfo_dir,
              args.year_table,
              args.source_strain)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
