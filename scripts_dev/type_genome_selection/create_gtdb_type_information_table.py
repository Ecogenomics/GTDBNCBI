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
import copy
from itertools import islice
from datetime import datetime
import math


class InfoGenerator(object):
    """Main class
      """

    def __init__(self, metadata_file, ncbi_names_file, lpsn_dir, dsmz_dir, straininfo_dir, output_dir, cpus=None):
        """Initialization."""
        self.metadata_dictionary, list_ids = self.load_metadata_dictionary(
            metadata_file)
        self.lpsn_strains_dic = self.load_lpsn_strains_dictionary(lpsn_dir)
        self.dsmz_strains_dic = self.load_dsmz_strains_dictionary(dsmz_dir)
        self.straininfo_strains_dic = self.load_straininfo_strains_dictionary(
            straininfo_dir)
        self.ncbi_names_dic = self.load_ncbi_names(ncbi_names_file, list_ids)
        self.output_dir = output_dir
        self.threads = int(cpus)
        self.pool_size = int(cpus)  # your "parallelness"
        self.pool = Pool(self.pool_size)

        self.startTime = datetime.now()

    def take(self, n, iterable):
        "Return first n items of the iterable as a list"
        return list(islice(iterable, n))

    def load_ncbi_names(self, ncbi_names_file, ids_of_interest):
        """ Load the Names.dmp file"""
        ncbi_names_dic = {}
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
                count += 1
        return ncbi_names_dic

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
                if not infos[gtdb_accession_index].startswith('U_'):
                    metadata_dictionary[infos[gtdb_accession_index]] = {
                        'ncbi_organism_name': infos[gtdb_ncbi_organism_name_index],
                        'taxonomy_species_name': infos[gtdb_taxonomy_species_name_index].split(';')[6].replace('s__',
                                                                                                               ''),
                        'strain_identifiers': infos[gtdb_strain_identifiers_index],
                        'ncbi_type_material_designation': infos[gtdb_ncbi_type_material_designation_index],
                        'ncbi_taxonomy_unfiltered': infos[gtdb_ncbi_taxonomy_unfiltered_index],
                        'ncbi_taxid': int(infos[gtdb_ncbi_taxid_index])}
                    list_taxid.append(int(infos[gtdb_ncbi_taxid_index]))
        return (metadata_dictionary, list(set(list_taxid)))

    def load_straininfo_strains_dictionary(self, straininfo_dir):
        # We load the dictionary of strains from Straininfo
        straininfo_strains_dic = {}
        p = re.compile('(\w+)\s\(now\s(\w+)\)\s(\d+)', re.IGNORECASE)
        with open(os.path.join(straininfo_dir, 'straininfo_strains.tsv')) as ststr:
            ststr.readline()
            for line in ststr:
                infos = line.rstrip('\n').split('\t')
                if len(infos) < 2:
                    print infos
                else:
                    new_ls = []
                    list_strains = infos[1].split()
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

                        straininfo_strains_dic[infos[0]] = new_ls
        return straininfo_strains_dic

    def load_dsmz_strains_dictionary(self, dsmz_dir):
        # We load the dictionary of strains from DSMZ
        dsmz_strains_dic = {}
        with open(os.path.join(dsmz_dir, 'dsmz_strains.tsv')) as dsstr:
            dsstr.readline()
            for line in dsstr:
                infos = line.rstrip('\n').split('\t')
                if len(infos) < 2:
                    print infos
                else:
                    dsmz_strains_dic[infos[0]] = infos[1]
        return dsmz_strains_dic

    def load_lpsn_strains_dictionary(self, lpsn_dir):
        # We load the dictionary of strains from LPSN
        lpsn_strains_dic = {}
        with open(os.path.join(lpsn_dir, 'lpsn_strains.tsv')) as lpstr:
            lpstr.readline()
            for line in lpstr:
                infos = line.rstrip('\n').split('\t')
                if len(infos) < 2:
                    print infos
                else:
                    lpsn_strains_dic[infos[0]] = infos[1]
        return lpsn_strains_dic

    def worker(self, mini_dict, i, out_q, strain_dictionary):
        count = 0
        # For each
        for acc, info_genomes in mini_dict.iteritems():
            count += 1
            # print "Thread {}: {}/{} done. ({}%)".format(i, count,
            # num_element, round((float(count) * 100 / num_element), 3))
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
            # if '[' and ']' starts the potentual name (i.e. st__[Eubacterium]
            # siraeum 70/3 ) those brackets are removed

            potential_names = [ncbi_organism_name]

            if ncbi_st_name:
                ncbi_taxonomy_name = ncbi_st_name
                if ncbi_taxonomy_name.startswith('['):
                    ncbi_taxonomy_name = ncbi_taxonomy_name.replace(
                        '[', '', 1).replace(']', '', 1)
                potential_names.append(ncbi_taxonomy_name)
            # if the last element is a subspecies
            if 'subsp.' in last_element:
                ncbi_taxonomy_name = last_element
                if ncbi_taxonomy_name.startswith('['):
                    ncbi_taxonomy_name = ncbi_taxonomy_name.replace(
                        '[', '', 1).replace(']', '', 1)
                potential_names.append(ncbi_taxonomy_name)
            # If the ncbi species name in the ncbi organism name we do not store it
            # otherwise ( like for RS_GCF_000469465.1 ) ncbi species name is
            # added to the potential list of names

            if ncbi_spe_name:
                p = re.compile(ncbi_spe_name + "(\s|\n)", re.IGNORECASE)
                matches = p.search(ncbi_organism_name)
                if matches is None:
                    ncbi_taxonomy_name = ncbi_spe_name
                    # if the ncbi organism name has a subspcies suffix, the suffix is added to the ncbi species name
                    # i.e ncbi organism name = Clavibacter michiganense subsp. insidiosum
                    # ncbi species name = Clavibacteria michiganense
                    # new potential name = Clavibacteria michiganense subsp.
                    # insidiosum
                    p = re.compile('subsp\. \w+', re.IGNORECASE)
                    matches = p.search(ncbi_organism_name)
                    if matches:
                        ncbi_taxonomy_name = "{} {}".format(
                            ncbi_taxonomy_name, matches.group(0))
                    if ncbi_taxonomy_name.startswith('['):
                        ncbi_taxonomy_name = ncbi_taxonomy_name.replace(
                            '[', '', 1).replace(']', '', 1)
                    potential_names.append(ncbi_taxonomy_name)
            # if the last element is binomial
            if not last_element.startswith('st__') and not last_element.startswith('s__') and len(last_element.split()) == 2:
                ncbi_taxonomy_name = last_element
                # if the ncbi organism name has a subspcies suffix, the suffix is added to the last_element
                # i.e ncbi organism name = Clavibacter michiganense subsp. insidiosum
                # last_element = Clavibacteria michiganense
                # new potential name = Clavibacteria michiganense subsp.
                # insidiosum
                p = re.compile('subsp\. \w+', re.IGNORECASE)
                matches = p.search(ncbi_organism_name)
                if matches:
                    ncbi_taxonomy_name = "{} {}".format(
                        ncbi_taxonomy_name, matches.group(0))
                if ncbi_taxonomy_name.startswith('['):
                    ncbi_taxonomy_name = ncbi_taxonomy_name.replace(
                        '[', '', 1).replace(']', '', 1)
                potential_names.append(ncbi_taxonomy_name)

                # we add all potential names from names.dmp
                if info_genomes.get('ncbi_taxid') in self.ncbi_names_dic:
                    potential_names.extend(self.ncbi_names_dic.get(
                        info_genomes.get('ncbi_taxid')))

            # We remove duplicates
            set_potential_names = set(potential_names)
            list_spes = []
            for pot_name in set_potential_names:
                # we look if any names in the lpsn strain dictionary is present in
                # the potential name
                for spe in strain_dictionary:
                    p = re.compile(spe + '(\s|\n)', re.IGNORECASE)
                    if 'subsp.' not in pot_name and p.search(pot_name):
                        list_spes.append(spe)
                        spe_list = spe.split()
                        # we recreate the subspecies name for a specie name
                        if len(spe_list) == 2:
                            subspe = "{} {} subsp. {}".format(
                                spe_list[0], spe_list[1], spe_list[1])
                            if subspe in strain_dictionary:
                                list_spes.append(subspe)
                    elif 'subsp.' in spe and p.search(pot_name):
                        list_spes.append(spe)
                        spe_list = spe.split()
                        # if the species name in the strain dictionary
                        # like 'Clavibacter michiganense subsp. michiganense'
                        # and 'Clavibacter michiganense' is also in the strain
                        # dictioanry we add both
                        if len(spe_list) == 4 and spe_list[1] == spe_list[3] \
                                and spe_list[0] + " " + spe_list[1] in strain_dictionary:
                            list_spes.append(spe_list[0] + " " + spe_list[1])

            istype = False

            # list_spes is the list of species name in LPSN,DSMZ or Straininfo
            # that have a match in the list of potential names

            if list_spes:
                set_spe = set(list_spes)
                for spe_name in set_spe:
                    list_strains = strain_dictionary.get(spe_name)
                    for strain in list_strains.split("="):
                        # "official" strains should have the formats "AAA123" or AAA 123"
                        p = re.compile('\w+\s?\d+')
                        if p.match(strain):
                            # if the strain if found in the list of potential
                            # names or strains_identifiers from NCBI, this
                            # genome is a type strain
                            if self.metadata_dictionary.get('strain_identifiers') is not None and (strain in self.metadata_dictionary.get('strain_identifiers') or strain in " ".join(set_potential_names)):
                                istype = True
                            elif strain in " ".join(set_potential_names):
                                istype = True
            out_q.put((acc, istype))
        return True

    def splitchunks(self, d, n):
        chunksize = int(math.ceil(len(d) / float(n)))
        it = iter(d)
        for _ in xrange(0, len(d), chunksize):
            yield {k: d[k] for k in islice(it, chunksize)}

    def parse_strains(self, strain_dictionary, filename):
        genome_strain_dict = {}
        manager = multiprocessing.Manager()
        out_q = manager.Queue()
        procs = []
        nprocs = self.threads
        # we split the metatadata dictionary to n subdictionaries of same length
        # n = number of processors
        for i, item in enumerate(self.splitchunks(self.metadata_dictionary, nprocs)):
            p = multiprocessing.Process(
                target=self.worker,
                args=(item, i, out_q, strain_dictionary))
            procs.append(p)
            p.start()

        # Collect all results into a single result dict. We know how many dicts
        # with results to expect.
        while out_q.empty():
            time.sleep(1)

        # Wait for all worker processes to finish
        results = {}

        for i in range(len(self.metadata_dictionary)):
            id, type = out_q.get()
            results[id] = type

        file_out = open(filename, 'w')

        for k, infos in results.iteritems():
            file_out.write("{}\t{}\n".format(k, infos))
        file_out.close()

    def run(self, sourcest):
        if sourcest == 'lpsn':
            self.parse_strains(self.lpsn_strains_dic, os.path.join(
                self.output_dir, 'lpsn_summary.txt'))
        elif sourcest == 'dsmz':
            self.parse_strains(self.dsmz_strains_dic, os.path.join(
                self.output_dir, 'dsmz_summary.txt'))
        elif sourcest == 'straininfo':
            self.parse_strains(self.lpsn_strains_dic, os.path.join(
                self.output_dir, 'straininfo_summary.txt'))
        elif sourcest == 'all':
            self.parse_strains(self.lpsn_strains_dic, os.path.join(
                self.output_dir, 'lpsn_summary.txt'))
            self.parse_strains(self.dsmz_strains_dic, os.path.join(
                self.output_dir, 'dsmz_summary.txt'))
            self.parse_strains(self.lpsn_strains_dic, os.path.join(
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
    parser.add_argument('--metadata_file',
                        help='Metadata file generated by gtdb')
    parser.add_argument('--ncbi_names', help='Names.dmp file from NCBI')
    parser.add_argument('--cpus',
                        help='Number of threads.')
    parser.add_argument('--output_dir',
                        help='Output directory.')
    parser.add_argument('--source_strain', choices=['all', 'lpsn', 'dsmz',
                                                    'straininfo'], default='all', help='select LPSN,Straininfo,DSMZ to parse')

    args = parser.parse_args()

    try:
        typeinfogenerator = InfoGenerator(
            args.metadata_file, args.ncbi_names, args.lpsn_dir, args.dsmz_dir, args.straininfo_dir, args.output_dir, args.cpus)
        typeinfogenerator.run(args.source_strain)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
