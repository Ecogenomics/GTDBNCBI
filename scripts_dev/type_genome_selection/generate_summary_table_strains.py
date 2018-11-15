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

__prog_name__ = 'generate_summary_table_strains.py'
__prog_desc__ = 'Generate a tab delimited file listing all information used to select type strains.'

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
import pickle
import os
from distutils.util import strtobool


class SummaryEditor(object):
    """Main class
      """

    def __init__(self, lpsn_strain_summary, dsmz_strain_summary, straininfo_strain_summary, metadata_file,
                 ncbi_names, ncbi_pickle, lpsn_type_genus_file, dsmz_type_genus_file):
        """Initialization."""
        self.metadata_dictionary, list_ids = self.load_metadata_dictionary(
            metadata_file)

        self.lpsn_strain_dict = self.load_strain_dic(lpsn_strain_summary)
        self.dsmz_strain_dict = self.load_strain_dic(dsmz_strain_summary)
        self.straininfo_strain_dict = self.load_strain_dic(
            straininfo_strain_summary)

        #print self.dsmz_strain_dict
        self.lpsn_type_genus_list = self.load_type_genus_file(lpsn_type_genus_file)
        self.dsmz_type_genus_list = self.load_type_genus_file(dsmz_type_genus_file)
        self.combine_list = list(self.dsmz_type_genus_list)
        self.combine_list.extend(self.lpsn_type_genus_list)
        self.ncbi_names_dic = {}
        # if the pickle dictionary exists ,we load it
        # if not it will be generated for the next time we run the script
        if os.path.isfile(ncbi_pickle):
            with open(ncbi_pickle, 'rb') as f:
                self.ncbi_names_dic = pickle.load(f)
        else:
            self.ncbi_names_dic = self.load_ncbi_names(ncbi_names, list_ids)
            with open(ncbi_pickle, 'wb') as f:
                pickle.dump(self.ncbi_names_dic, f, pickle.HIGHEST_PROTOCOL)

    def load_type_genus_file(self,genus_file):
        list_type_genus_species = []
        with open(genus_file) as gf:
            headers = gf.readline()
            for line in gf:
                infos = line.rstrip('\n').split('\t')
                if infos[1] != '':
                    list_type_genus_species.append(infos[0].replace('s__',''))
        return list_type_genus_species

    def load_strain_dic(self, strain_summary_file):
        """
        LPSN,DSMZ,Straininfo are all the same format
        (genome_id \t true/false \t true/false)
        """
        result_dict = {}
        with open(strain_summary_file) as ssf:
            for line in ssf:
                infos = line.rstrip().split('\t')
                if len(infos) == 2 :
                    result_dict[infos[0]] = {"type":strtobool(infos[1]),'os':strtobool(infos[3]),'neotype':False}
                else:
                    result_dict[infos[0]] = {"type":strtobool(infos[1]),'os':strtobool(infos[3]),'neotype':strtobool(infos[2])}
        return result_dict

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

    def load_ncbi_names(self, ncbi_names_file, ids_of_interest):
        """ Parses the dmp file """
        ncbi_names_dic = {}
        print "Loading NCBI names."
        num_lines = sum(1 for line in open(ncbi_names_file))
        count = 0
        with open(ncbi_names_file) as nnf:
            for line in nnf:
                count += 1
                sys.stdout.write("{}% complete\r".format(
                    round((float(count) * 100 / num_lines), 3)))
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
        return ncbi_names_dic

    def run(self, outfile):
        outf = open(outfile, 'w')
        out_syn = open(os.path.join(os.path.dirname(outfile),'only_synonyms.tsv'),'w')

        outf.write(
            "accession\tgtdb_type_material\tgtdb_type_material_sources\t")
        outf.write(
            "gtdb_type_species_of_genus\tgtdb_type_strain_of_species\tgtdb_type_strain_of_subspecies\t")
        outf.write('ncbi_names\n')

        count = 0
        for acc, infos_genomes in self.metadata_dictionary.iteritems():
            infores = []
            infores.append(acc)
            sources_type = []
            if self.lpsn_strain_dict.get(acc).get('neotype'):
                sources_type.append('lpsn(neotype)')
            elif self.lpsn_strain_dict.get(acc).get('type'):
                sources_type.append('lpsn')
            #print self.dsmz_type_genus_list
            if self.dsmz_strain_dict.get(acc).get('type'):
                sources_type.append('dsmz')
            if self.straininfo_strain_dict.get(acc).get('type'):
                sources_type.append('straininfo')
            #if infos_genomes.get('ncbi_type_material_designation') is not None and infos_genomes.get('ncbi_type_material_designation') != '' and infos_genomes.get('ncbi_type_material_designation') != 'none':
            #    sources_type.append('ncbi')
            if len(sources_type) > 0 :
                infores.append('true')
            else:
                infores.append('false')
            infores.append('_'.join(sources_type))
            #TODO DO LATER : Implement type material of a genus
            infores.append('false')
            ######
            type_strain = False
            if len(sources_type) > 0 :
                if 'subsp.' not in infos_genomes.get('ncbi_organism_name'):
                    infores.append('true')
                    type_strain = True

                else:
                    spe_list = infos_genomes.get('ncbi_organism_name').split()
                    subsp_index = spe_list.index('subsp.')
                    if spe_list[subsp_index -1] == spe_list[subsp_index +1]:
                        infores.append('true')
                        type_strain = True
                    else:
                        infores.append('false')
            else:
                infores.append('false')
            if type_strain and any(spename in infos_genomes.get('ncbi_organism_name')
                                   for spename in self.combine_list):
                infores[3] = 'true'
            if len(sources_type) > 0 and 'subsp.' in infos_genomes.get('ncbi_organism_name'):
                spe_list = infos_genomes.get('ncbi_organism_name').split()
                subsp_index = spe_list.index('subsp.')
                if spe_list[subsp_index - 1] != spe_list[subsp_index + 1]:
                    infores.append('true')
                else:
                    infores.append('false')
            else:
                infores.append('false')
            #infores.append(infos_genomes.get('ncbi_organism_name'))
            if infos_genomes.get('ncbi_taxid') in self.ncbi_names_dic:
                infores.append(
                    '/'.join(self.ncbi_names_dic.get(infos_genomes.get('ncbi_taxid'))))
            else:
                infores.append("None")

            outf.write("{}\n".format("\t".join(infores)))

            if 'true' in infores and self.lpsn_strain_dict.get(acc).get('os') and self.dsmz_strain_dict.get(acc).get('os') and self.straininfo_strain_dict.get(acc).get('os'):
                print infos_genomes.keys()
                out_syn.write('{}\t{}\t{}\t{}\n'.format(acc,infos_genomes.get('ncbi_organism_name'),
                                                        infos_genomes.get('strain_identifiers'),
                              '/'.join(self.ncbi_names_dic.get(infos_genomes.get('ncbi_taxid')))))


        outf.close()


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--lpsn_strain_summary', help='LPSN strain summary file created by the create_gtdb_type_information_table.py script.')
    parser.add_argument(
        '--dsmz_strain_summary', help='DSMZ strain summary file created by the create_gtdb_type_information_table.py script.')
    parser.add_argument(
        '--straininfo_strain_summary', help='Straininfo strain summary file created by the create_gtdb_type_information_table.py script.')
    parser.add_argument('--lpsn_type_genus_file', help='LPSN genus summary file created by previous step')
    parser.add_argument('--dsmz_type_genus_file', help='DSMZ genus summary file created by previous step')
    parser.add_argument('--metadata_file',
                        help='Metadata file generated by gtdb')
    parser.add_argument('--ncbi_names', help='Names.dmp file from NCBI')
    parser.add_argument('--ncbi_pickle',
                        help='prepocessed dictionary Names.dmp file from NCBI. If the file doesnt exist , It will be created for future run of the script')
    parser.add_argument('--output',
                        help='Output file.')

    args = parser.parse_args()

    try:
        summaryeditor = SummaryEditor(
            args.lpsn_strain_summary, args.dsmz_strain_summary,
            args.straininfo_strain_summary, args.metadata_file,
            args.ncbi_names, args.ncbi_pickle, args.lpsn_type_genus_file, args.dsmz_type_genus_file)
        summaryeditor.run(args.output)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
