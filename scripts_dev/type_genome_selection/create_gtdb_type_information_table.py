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

__prog_name__ = 'create_gtdb_type_information_table.py'
__prog_desc__ = 'Parse multiple sources to pick type species,type subspecies and type genus'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2018'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

import sys
import argparse
import os
import re
from multiprocessing.pool import ThreadPool as Pool
import multiprocessing
import time
from itertools import islice
from datetime import datetime
import math
import pickle
from collections import Counter
from fuzzywuzzy import fuzz
from fuzzywuzzy import process


class InfoGenerator(object):
    """Main class
      """

    def __init__(self, metadata_file, ncbi_names_file, lpsn_dir, dsmz_dir, straininfo_dir, ncbi_pickle, output_dir, year_table, cpus=None):
        """Initialization."""
        self.metadata_dictionary, list_ids = self.load_metadata_dictionary(
            metadata_file)
        self.year_table = self.load_year_dict(year_table)
        self.lpsn_strains_dic = self.load_lpsn_strains_dictionary(lpsn_dir)
        self.dsmz_strains_dic = self.load_dsmz_strains_dictionary(dsmz_dir)
        self.straininfo_strains_dic = self.load_straininfo_strains_dictionary(
            straininfo_dir)
        self.ncbi_names_dic = {}
        # if the pickle dictionary exists ,we load it
        # if not it will be generated for the next time we run the script
        if os.path.isfile(ncbi_pickle):
            with open(ncbi_pickle, 'rb') as f:
                self.ncbi_names_dic = pickle.load(f)
        else:
            self.ncbi_names_dic = self.load_ncbi_names(
                ncbi_names_file, list_ids)
            with open(ncbi_pickle, 'wb') as f:
                pickle.dump(self.ncbi_names_dic, f, pickle.HIGHEST_PROTOCOL)

        for k, v in self.ncbi_names_dic.iteritems():
            if len(set(v.get('synonym')).intersection(v.get('scientific name'))) > 0 or len(set(v.get('synonym')).intersection(v.get('equivalent name'))) > 0:
                print 'ERROR'
                print v.get('synonym')
                print v.get('scientific name')
                print v.get('equivalent name')

        self.output_dir = output_dir
        self.threads = int(cpus)
        self.pool_size = int(cpus)  # your "parallelness"
        self.pool = Pool(self.pool_size)

        self.startTime = datetime.now()

    def take(self, n, iterable):
        "Return first n items of the iterable as a list"
        return list(islice(iterable, n))

    def load_year_dict(self, year_table):
        dict_date = {}
        with open(year_table) as yt:
            for line in yt:
                infos = line.rstrip('\n').split('\t')
                dict_date[infos[0]] = {'lpsn': infos[1],
                                       'dsmz': infos[2], 'straininfo': infos[3]}
        return dict_date

    def load_ncbi_names(self, ncbi_names_file, ids_of_interest):
        """ Load the Names.dmp file"""
        ncbi_names_dic = {}
        category_names = {}
        print "Loading NCBI names."
        num_lines = sum(1 for line in open(ncbi_names_file))
        count = 1
        with open(ncbi_names_file) as nnf:
            for line in nnf:
                sys.stdout.write("{}% complete\r".format(
                    round((float(count) * 100 / num_lines), 3)))
                # sys.stdout.write("{}/{}\r".format(count, num_lines))
                sys.stdout.flush()
                infos = line.split('|')
                current_id = int(infos[0].strip())
                if current_id in ids_of_interest:
                    if infos[3].strip() in ['misspelling', 'synonym', 'equivalent name', 'scientific name']:
                        if current_id in ncbi_names_dic:
                            ncbi_names_dic.get(current_id).append(
                                infos[1].strip())
                        else:
                            ncbi_names_dic[current_id] = [infos[1].strip()]
                        if current_id not in category_names:
                            category_names[current_id] = {'misspelling': [], 'synonym': [
                            ], 'equivalent name': [], 'scientific name': []}

                        category_names.get(current_id).get(
                            infos[3].strip()).append(infos[1].strip())
                count += 1
        return category_names

    def load_metadata_dictionary(self, metadata_file):
        """
        We are parsing as much data as possible from the metadata file
        """
        metadata_dictionary = {}
        with open(metadata_file) as metaf:
            headers_line = metaf.readline()
            if '\t' in headers_line:
                separator_info = '\t'
            else:
                separator_info = ','
            headers = headers_line.rstrip('\n').split(separator_info)

            gtdb_ncbi_organism_name_index = headers.index('ncbi_organism_name')
            gtdb_ncbi_type_material_designation_index = headers.index(
                'ncbi_type_material_designation')
            gtdb_accession_index = headers.index('accession')
            gtdb_taxonomy_species_name_index = headers.index('ncbi_taxonomy')
            gtdb_strain_identifiers_index = headers.index(
                'ncbi_strain_identifiers')
            gtdb_ncbi_taxonomy_unfiltered_index = headers.index(
                'ncbi_taxonomy_unfiltered')
            gtdb_ncbi_taxid_index = headers.index(
                'ncbi_taxid')

            list_taxid = []

            for line in metaf:
                infos = line.rstrip('\n').split(separator_info)
                # and infos[gtdb_accession_index] == 'RS_GCF_001975225.1':
                if not infos[gtdb_accession_index].startswith('U_'):
                    # We processed the strains to be 'standardised'
                    pattern = re.compile('[\W_]+')
                    created_list = re.split(
                        ';|/', infos[gtdb_strain_identifiers_index])

                    list_strains = [pattern.sub('', a.strip()).upper(
                    ) for a in created_list if (a != '' and a != 'none')]
                    metadata_dictionary[infos[gtdb_accession_index]] = {
                        'ncbi_organism_name': infos[gtdb_ncbi_organism_name_index],
                        'taxonomy_species_name': infos[gtdb_taxonomy_species_name_index].split(';')[6].replace('s__',
                                                                                                               ''),
                        'strain_identifiers': set(list_strains),
                        'ncbi_type_material_designation': infos[gtdb_ncbi_type_material_designation_index],
                        'ncbi_taxonomy_unfiltered': infos[gtdb_ncbi_taxonomy_unfiltered_index],
                        'ncbi_taxid': int(infos[gtdb_ncbi_taxid_index])}
                    list_taxid.append(int(infos[gtdb_ncbi_taxid_index]))
        return (metadata_dictionary, list(set(list_taxid)))

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

    def remove_brackets(self, name):
        if name.startswith('['):
            name = name.replace('[', '', 1).replace(']', '', 1)
        return name

    def add_taxonomic_names(self, potential_names, ncbi_organism_name, ncbi_spe_name, ncbi_st_name, last_element):
        if ncbi_st_name:
            ncbi_taxonomy_name = self.remove_brackets(ncbi_st_name)
            potential_names.append(ncbi_taxonomy_name)
        # if the last element is a subspecies
        if 'subsp.' in last_element:
            ncbi_taxonomy_name = self.remove_brackets(last_element)
            potential_names.append(ncbi_taxonomy_name)
        # If the ncbi species name in the ncbi organism name we do not store it
        # otherwise ( like for RS_GCF_000469465.1 ) ncbi species name is
        # added to the potential list of names
        if ncbi_spe_name:
            p = re.compile(ncbi_spe_name + "(\s|\n)", re.IGNORECASE)
            matches = p.search(ncbi_organism_name)
            if matches is None:
                ncbi_taxonomy_name = ncbi_spe_name
                # if the ncbi organism name has a subspecies suffix, the suffix is added to the ncbi species name
                # i.e ncbi organism name = Clavibacter michiganense subsp. insidiosum
                # ncbi species name = Clavibacteria michiganense
                # new potential name = Clavibacteria michiganense subsp.
                # insidiosum
                p = re.compile('subsp\. \w+', re.IGNORECASE)
                matches = p.search(ncbi_organism_name)
                if matches:
                    ncbi_taxonomy_name = "{} {}".format(
                        ncbi_taxonomy_name, matches.group(0))
                ncbi_taxonomy_name = self.remove_brackets(
                    ncbi_taxonomy_name)
                potential_names.append(ncbi_taxonomy_name)
        # if the last element is binomial
        if not last_element.startswith('st__') and not last_element.startswith('s__') and len(last_element.split()) == 2:
            ncbi_taxonomy_name = last_element
            # if the ncbi organism name has a subspecies suffix, the suffix is added to the last_element
            # i.e ncbi organism name = Clavibacter michiganense subsp. insidiosum
            # last_element = Clavibacteria michiganense
            # new potential name = Clavibacteria michiganense subsp.
            # insidiosum
            p = re.compile('subsp\. \w+', re.IGNORECASE)
            matches = p.search(ncbi_organism_name)
            if matches:
                ncbi_taxonomy_name = "{} {}".format(
                    ncbi_taxonomy_name, matches.group(0))
            ncbi_taxonomy_name = self.remove_brackets(ncbi_taxonomy_name)
            potential_names.append(ncbi_taxonomy_name)
        return potential_names

    def worker(self, mini_dict, sourcest, i, out_q, strain_dictionary):
        #count = 0
        #num_element = len(mini_dict)
        start = time.time()
        for acc, info_genomes in mini_dict.iteritems():
            #==================================================================
            # if count > 100:
            #    break
            # count += 1
            # print "Thread {}: {}/{} done. ({}%)".format(i, count,
            #                                             num_element, round((float(count) * 100 / num_element), 3))
            #==================================================================
            ncbi_organism_name = info_genomes.get('ncbi_organism_name')
            ncbi_unfiltered_tax = info_genomes.get('ncbi_taxonomy_unfiltered')
            ncbi_unfiltered_list = ncbi_unfiltered_tax.split(';')
            ncbi_taxonomy_name = None
            last_element = ncbi_unfiltered_list[-1].split('__')[1]
            ncbi_spe_name = None
            ncbi_st_name = None

            for item in ncbi_unfiltered_list:
                if item.startswith('s__'):
                    ncbi_spe_name = item[3:]
                elif item.startswith('st__'):
                    ncbi_st_name = item[4:]
                    break

            # we store the last_element, the species name and the strain name from the
            # ncbi_unfiltered taxonomy field
            # if '[' and ']' starts the potential name (i.e. st__[Eubacterium]
            # siraeum 70/3 ) those brackets are removed

            potential_names = [ncbi_organism_name]
            potential_names = self.add_taxonomic_names(
                potential_names, ncbi_organism_name, ncbi_spe_name, ncbi_st_name, last_element)

            # scientific names are considered official names
            if info_genomes.get('ncbi_taxid') in self.ncbi_names_dic:
                potential_names.extend(self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('scientific name'))
            official_potential_names = set(potential_names)

            # we add all potential names from names.dmp
            #['misspelling', 'synonym', 'equivalent name']
            misspelling_names = []
            synonyms = []
            equivalent_names = []
            scientific_names = []
            potential_names_notoff = []

            if info_genomes.get('ncbi_taxid') in self.ncbi_names_dic:
                potential_names_notoff.extend(self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('misspelling'))
                potential_names_notoff.extend(self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('synonym'))
                potential_names_notoff.extend(self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('equivalent name'))

                misspelling_names = self.standardise_names(self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('misspelling'), strain_dictionary)
                unprocessed_misspelling_names = self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('misspelling')
                synonyms = self.standardise_names(self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('synonym'), strain_dictionary)
                unprocessed_synonyms = self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('synonym')
                equivalent_names = self.standardise_names(self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('equivalent name'), strain_dictionary)
                unprocessed_equivalent_names = self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('equivalent name')
                scientific_names = self.standardise_names(self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('scientific name'), strain_dictionary)
                unprocessed_scientific_names = self.ncbi_names_dic.get(
                    info_genomes.get('ncbi_taxid')).get('scientific names')

            # We remove duplicates
            non_official_potential_names = set(potential_names_notoff)

            official_spe_names = self.standardise_names(
                official_potential_names, strain_dictionary)
            non_official_spe_names = self.standardise_names(
                non_official_potential_names, strain_dictionary)

            category_information_list = self.strain_match(
                acc, official_spe_names, official_potential_names, misspelling_names, synonyms, equivalent_names, strain_dictionary, sourcest, True)
            if len(category_information_list) == 0 or category_information_list[0][0] != 'official_name':
                category_information_list = self.strain_match(
                    acc, non_official_spe_names, non_official_potential_names, misspelling_names, synonyms, equivalent_names, strain_dictionary, sourcest, False)

            end = time.time()
            # print "Mid:{}\tLast:{}".format(mid - start, end - mid)

            for item_list in category_information_list:
                category_name = item_list[0]
                istype = item_list[1]
                isneotype = item_list[2]
                spename_off = item_list[3]
                year_date = item_list[4]

                # dereplication of official names
                dereplicated_names_list = []
                fuzzy_score = 0
                fuzzy_match = ''
                for name in official_spe_names:
                    if len(name.split(' ')) == 4:
                        infos = name.split(' ')
                        if infos[1] != infos[3]:
                            dereplicated_names_list.append(name)
                            if fuzz.ratio(name, spename_off) > fuzzy_score:
                                fuzzy_score = fuzz.ratio(name, spename_off)
                                fuzzy_match = '{}/{}'.format(name, spename_off)
                        else:
                            dereplicated_names_list.append(
                                infos[0] + ' ' + infos[1])
                            if fuzz.ratio(infos[0] + ' ' + infos[1], spename_off) > fuzzy_score:
                                fuzzy_score = fuzz.ratio(
                                    infos[0] + ' ' + infos[1], spename_off)
                                fuzzy_match = '{}/{}'.format(
                                    infos[0] + ' ' + infos[1], spename_off)
                    else:
                        dereplicated_names_list.append(name)
                        if fuzz.ratio(name, spename_off) > fuzzy_score:
                            fuzzy_score = fuzz.ratio(name, spename_off)
                            fuzzy_match = '{}/{}'.format(name, spename_off)
                dereplicated_names = set(dereplicated_names_list)
                out_q.put((acc, year_date, istype, isneotype, category_name, spename_off, dereplicated_names, fuzzy_match, fuzzy_score, unprocessed_misspelling_names,
                           unprocessed_synonyms, unprocessed_equivalent_names, unprocessed_scientific_names))
                # out_q.put(None)
            start = time.time()
        return True

    def strains_iterate(self, acc, spe_name, list_strains, raw_potential_names, misspelling_names, synonyms, equivalent_names, isofficial, sourcest):
        istype = False
        year_date = ''
        category_name = ''
        spename_off = ''
        for strain in list_strains.split("="):
            if len(strain) > 1:
                # if the strain if found in the list of potential
                # names or strains_identifiers from NCBI, this
                # genome is a type strain
                pattern = re.compile('[\W_]+')
                standard_pot_names = {pattern.sub(
                    '', a).upper(): a for a in raw_potential_names}
                if self.metadata_dictionary.get(acc).get('strain_identifiers') is not None and strain in self.metadata_dictionary.get(acc).get('strain_identifiers'):
                    istype = True
                    spename_off = spe_name
                    if not isofficial:
                        category_name = self.select_category_name(
                            spe_name, misspelling_names, synonyms, equivalent_names)
                    else:
                        if spe_name in self.year_table:
                            if self.year_table.get(spe_name).get(sourcest) != '':
                                if (year_date != '' and int(year_date) < int(self.year_table.get(spe_name).get(sourcest))) or year_date == '':
                                    year_date = self.year_table.get(
                                        spe_name).get(sourcest)
                        category_name = 'official_name'

                else:
                    for standard_pot_name in standard_pot_names:
                        if strain in standard_pot_name:
                            first_char = strain[0]
                            last_char = strain[-1]
                            p = re.compile(' {}'.format(
                                first_char), re.IGNORECASE)
                            matches_beginning = p.search(
                                standard_pot_names.get(standard_pot_name))
                            q = re.compile('{}(\s|$)'.format(
                                last_char), re.IGNORECASE)
                            matches_end = q.search(
                                standard_pot_names.get(standard_pot_name))
                            if matches_beginning and matches_end:
                                if not isofficial:
                                        istype = True
                                        spename_off = spe_name
                                        category_name = self.select_category_name(
                                            spe_name, misspelling_names, synonyms, equivalent_names)
                                else:
                                    istype = True
                                    spename_off = spe_name
                                    if spe_name in self.year_table:
                                        if self.year_table.get(spe_name).get(sourcest) != '':
                                            if (year_date != '' and int(year_date) < int(self.year_table.get(spe_name).get(sourcest))) or year_date == '':
                                                year_date = self.year_table.get(
                                                    spe_name).get(sourcest)
                                    category_name = 'official_name'
        return (category_name, istype, spename_off, year_date)

    def strain_match(self, acc, list_names, raw_potential_names, misspelling_names, synonyms, equivalent_names, strain_dictionary, sourcest, isofficial):
        # list_names is the list of species name in LPSN,DSMZ or Straininfo
        # that have a match in the list of potential names

        list_category = []

        if list_names:
            # we remove duplicates a second time after standardisation of the
            # names
            set_spe = set(list_names)
            for spe_name in set_spe:
                istype = False
                isneotype = False
                year_date = ''
                category_name = ''
                spename_off = ''
                if spe_name in strain_dictionary:
                    list_strains = []
                    if sourcest == 'lpsn':
                        # lpsn as information for both strains and neotype
                        # strains
                        list_strains = strain_dictionary.get(
                            spe_name).get('strains')
                    else:
                        list_strains = strain_dictionary.get(spe_name)
                    category_name, istype, spename_off, year_date = self.strains_iterate(
                        acc, spe_name, list_strains, raw_potential_names, misspelling_names, synonyms, equivalent_names, isofficial, sourcest)
                    if sourcest == 'lpsn':
                        list_neotypes = strain_dictionary.get(
                            spe_name).get('neotypes')
                        _neotype_category_name, isneotype, _neotype_spename_off, _neotype_year_date = self.strains_iterate(
                            acc, spe_name, list_neotypes, raw_potential_names, misspelling_names, synonyms, equivalent_names, isofficial, sourcest)
                if category_name == 'official_name':
                    list_category=[[category_name, istype, isneotype, spename_off, year_date]]
                    return list_category
                elif category_name != '' :
                    list_category.append(
                        [category_name, istype, isneotype, spename_off, year_date])
                elif self.metadata_dictionary.get(acc).get('ncbi_type_material_designation')!='none':
                    list_category.append(
                        [category_name, istype, isneotype, spename_off, year_date])
        return list_category

    def select_category_name(self, spe_name, misspelling_names, synonyms, equivalent_names):
        if spe_name in misspelling_names:
            return 'misspelling name'
        elif spe_name in equivalent_names:
            return 'equivalent name'
        elif spe_name in synonyms:
            return 'synonyms'

    def standardise_names(self, potential_names, strain_dictionary):
        list_spes = []
        for pot_name in potential_names:
            pot_name_list = pot_name.split(' ')
            # if it's not a binomial name
            if len(pot_name_list) != 2:
                # we 'clean' the potential name and again if it's binomial
                pot_name = re.sub(r'(?i)(candidatus\s)', r'', pot_name)
                pot_name = re.sub(r'(?i)(serotype.*)', r'', pot_name)
                pot_name = re.sub(r'(?i)(ser\..*)', r'', pot_name)
                pot_name = re.sub(r'\"|\'', r'', pot_name)
                pot_name = self.remove_brackets(pot_name)
                pot_name = pot_name.strip()
                pot_name_list = pot_name.split(' ')
                if len(pot_name_list) == 4 and pot_name_list[2] == 'subsp.':
                    list_spes.append(pot_name)
                    # if the subspecies matches the species name
                    if pot_name_list[1] == pot_name_list[3]:
                        list_spes.append('{} {}'.format(
                            pot_name_list[0], pot_name_list[1]))
                # if the name is longer than 4 words but the 3rd word still
                # subsp, we assume than the 4 first words are a subspecies name
                elif len(pot_name_list) >= 4 and pot_name_list[2] == 'subsp.':
                    if pot_name_list[1] == pot_name_list[3]:
                        list_spes.append(
                            pot_name_list[0] + ' ' + pot_name_list[1])
                        list_spes.append('{0} {1} subsp. {1}'.format(
                            pot_name_list[0], pot_name_list[1]))
                    if pot_name_list[1] != pot_name_list[3]:
                        list_spes.append('{} {} {} {}'.format(
                            pot_name_list[0], pot_name_list[1], pot_name_list[2], pot_name_list[3]))
                elif len(pot_name_list) >= 2:
                    list_spes.append('{0} {1}'.format(
                        pot_name_list[0], pot_name_list[1]))
                    list_spes.append('{0} {1} subsp. {1}'.format(
                        pot_name_list[0], pot_name_list[1]))
            elif pot_name in strain_dictionary:
                list_spes.append(pot_name)
                list_spes.append('{0} {1} subsp. {1}'.format(
                    pot_name_list[0], pot_name_list[1]))
        return list_spes

    def splitchunks(self, d, n):
        chunksize = int(math.ceil(len(d) / float(n)))
        it = iter(d)
        for _ in xrange(0, len(d), chunksize):
            yield {k: d[k] for k in islice(it, chunksize)}

    def parse_strains(self, sourcest, strain_dictionary, filename):
        manager = multiprocessing.Manager()
        out_q = manager.Queue()
        procs = []
        nprocs = self.threads
        # we split the metatadata dictionary to n subdictionaries of same length
        # n = number of processors
        for i, item in enumerate(self.splitchunks(self.metadata_dictionary, nprocs)):
            p = multiprocessing.Process(
                target=self.worker,
                args=(item, sourcest, i, out_q, strain_dictionary))
            procs.append(p)
            p.start()

        # Collect all results into a single result dict. We know how many dicts
        # with results to expect.
        while out_q.empty():
            time.sleep(1)

        # Wait for all worker processes to finish
        results = {}
        # for i in range(len(self.metadata_dictionary)):
        count = 0
        isemptylist = False
        while not isemptylist:
            id_genome, year_date, type_strain, neotype, category_name, spename_off, dereplicated_names, fuzzy_match, fuzzy_score, misspelling_names, synonyms, equivalent_names, scientific_names = out_q.get()
            count += 1
            if id_genome in results:
                results.get(id_genome).append({
                    'type_strain': type_strain, 'yd': year_date, 'neotype': neotype, 'cat_name': category_name,
                    'species_name_match': spename_off, 'scientific_names': scientific_names, 'derep_set': dereplicated_names, 'fuzzy_match': fuzzy_match, 'fuzzy_score': fuzzy_score,
                    'synonyms': synonyms, 'misspellings': misspelling_names, 'equivalent_names': equivalent_names})
            else:
                results[id_genome] = [{
                    'type_strain': type_strain, 'yd': year_date, 'neotype': neotype, 'cat_name': category_name,
                    'species_name_match': spename_off, 'scientific_names': scientific_names, 'derep_set': dereplicated_names, 'fuzzy_match': fuzzy_match, 'fuzzy_score': fuzzy_score,
                    'synonyms': synonyms, 'misspellings': misspelling_names, 'equivalent_names': equivalent_names}]
            if out_q.empty():
                time.sleep(10)
                if out_q.empty():
                    isemptylist = True

        print "we process the data"
        file_catout = open(filename.replace('summary', 'category'), 'w')
        file_catout.write(
            'genome\tncbi_type_designation\tofficial_names\tmissspellings\tequivalent_names\tsynonyms\t{0}_match_type\t{0}_match_name\t{0}_fuzzy_score\t{0}_fuzzy_match\ttype_strain\tneotype\tpriority_date\n'.format(sourcest))
        for k, infos_list in results.iteritems():
            for infos in infos_list:
                if infos.get('cat_name') != '':
                    file_catout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        k,self.metadata_dictionary.get(k).get('ncbi_type_material_designation'), '/'.join(infos.get('derep_set')),
                        '/'.join(
                            infos.get('misspellings')),
                        '/'.join(
                            infos.get('equivalent_names')),
                        '/'.join(
                            infos.get('synonyms')),
                        infos.get(
                            'cat_name'),
                        infos.get(
                            'species_name_match'),
                        infos.get(
                            'fuzzy_score'),
                        infos.get(
                            'fuzzy_match'),
                        infos.get(
                            'type_strain'),
                        infos.get(
                            'neotype'),
                        infos.get('yd')))

        file_catout.close()

    def run(self, sourcest):
        if sourcest == 'lpsn':
            self.parse_strains(sourcest, self.lpsn_strains_dic, os.path.join(
                self.output_dir, 'lpsn_summary.txt'))
        elif sourcest == 'dsmz':
            self.parse_strains(sourcest, self.dsmz_strains_dic, os.path.join(
                self.output_dir, 'dsmz_summary.txt'))
        elif sourcest == 'straininfo':
            self.parse_strains(sourcest, self.straininfo_strains_dic, os.path.join(
                self.output_dir, 'straininfo_summary.txt'))
        elif sourcest == 'all':
            self.parse_strains('lpsn', self.lpsn_strains_dic, os.path.join(
                self.output_dir, 'lpsn_summary.txt'))
            self.parse_strains('dsmz', self.dsmz_strains_dic, os.path.join(
                self.output_dir, 'dsmz_summary.txt'))
            self.parse_strains('straininfo', self.straininfo_strains_dic, os.path.join(
                self.output_dir, 'straininfo_summary.txt'))


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--lpsn_dir', help='Directory including the 3 lpsn result files (lpsn_genera.tsv, lpsn_species.tsv and lpsn_strains.tsv ).')
    parser.add_argument(
        '--dsmz_dir', help='Directory including the 3 dsmz result files (dsmz_strains.tsv, dsmz_species.tsv and dsmz_genera.tsv ).')
    parser.add_argument(
        '--straininfo_dir', help='Directory including the straininfo result file (straininfo_strains.tsv).')
    parser.add_argument(
        '--year_table', help='Date table generated by generate_date_table.py .')
    parser.add_argument('--metadata_file',
                        help='Metadata file generated by gtdb')
    parser.add_argument('--ncbi_names', help='Names.dmp file from NCBI')
    parser.add_argument('--cpus',
                        help='Number of threads.')
    parser.add_argument('--output_dir',
                        help='Output directory.')
    parser.add_argument('--source_strain', choices=['all', 'lpsn', 'dsmz',
                                                    'straininfo'], default='all', help='select LPSN,Straininfo,DSMZ to parse')
    parser.add_argument('--ncbi_pickle',
                        help='prepocessed dictionary Names.dmp file from NCBI. If the file doesnt exist , It will be created for future run of the script')

    args = parser.parse_args()

    try:
        typeinfogenerator = InfoGenerator(
            args.metadata_file, args.ncbi_names, args.lpsn_dir, args.dsmz_dir, args.straininfo_dir, args.ncbi_pickle, args.output_dir, args.year_table, args.cpus)
        typeinfogenerator.run(args.source_strain)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise