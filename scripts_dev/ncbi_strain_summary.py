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

__prog_name__ = 'ncbi_strain_summary.py'
__prog_desc__ = 'Parse the assembly report file, the genomic.gbff file and the wgsmaster.gbff'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import os
import sys
import argparse
import tempfile
from collections import defaultdict
import re


class StrainParser(object):
    """Extract genes in nucleotide space."""

    def __init__(self, assembly_summary_genbank, assembly_summary_refseq):
        self.genbank_dictionary = self.parse_summary(assembly_summary_genbank)
        self.refseq_dictionary = self.parse_summary(assembly_summary_refseq)

    def file_len(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def parse_summary(self, assembly_summary):
        assembly_summary_dict = {}
        with open(assembly_summary) as as_file:
            as_file.readline()
            for line in as_file:
                if line.startswith('#'):
                    line = line.replace('# ', '')
                    headers = line.strip('\n').split('\t')
                    index_genome_id = headers.index('assembly_accession')
                    relation_to_type_material_index = headers.index(
                        'relation_to_type_material')
                else:
                    line_infos = line.strip('\n').split('\t')
                    assembly_summary_dict[line_infos[index_genome_id]
                                          ] = line_infos[relation_to_type_material_index]
        return assembly_summary_dict

    def run(self, genome_dir_file, out_dir):
        outf = open(os.path.join(out_dir, "strain_summary_file.txt"), "w")
        outf.write(
            "genome_id\tOrganism name\tstrain_identifiers\ttype_material_designation\n")
        pattern_strain = re.compile('^\s+\/strain=".+')
        pattern_isolate = re.compile('^\s+\/isolate=".+')
        number_of_genomes = self.file_len(genome_dir_file)
        count = 1
        with open(genome_dir_file, 'r') as genomelistfile:
            for line in genomelistfile:
                sys.stdout.write("{}% complete\r".format(
                    round((float(count) * 100 / number_of_genomes), 3)))
                # sys.stdout.write("{}/{}\r".format(count, num_lines))
                sys.stdout.flush()
                count += 1
                line_split = line.strip().split('\t')

                genome_id = line_split[0]
                #==============================================================
                # if genome_id not in ['GCF_000756285.1', 'GCF_900104035.1', 'GCF_000013925.1', 'GCF_900099745.1', 'GCF_000006605.1', 'GCF_000953375.1', 'GCF_000285595.1',
                #                      'GCF_001013905.1', 'GCF_001281105.1', 'GCF_000006805.1', 'GCF_000022905.1', 'GCF_002151505.1']:
                #     continue
                #==============================================================

                #==============================================================
                # if genome_id not in ['GCF_000006805.1']:
                #     continue
                #==============================================================

                genome_path = line_split[1]
                genome_dir_id = os.path.basename(os.path.normpath(genome_path))

                species = ''
                strain = ''
                isolate = ''
                typemat = ''
                wstrain = ''
                gstrains = []
                gstrain = ""
                wsisolate = ''
                gisolates = []
                gisolate = ""

                genbank_file = os.path.join(
                    line_split[1], genome_dir_id + '_assembly_report.txt')
                if os.path.exists(genbank_file):
                    with open(genbank_file, 'r') as gfile:
                        for gline in gfile:
                            if gline.startswith('# Organism name: '):
                                speline = gline.replace("# Organism name:", "")
                                species = re.sub(
                                    r'\([^)]*\)', '', speline).strip()
                            if gline.startswith('# Infraspecific name:'):
                                strline = gline.replace(
                                    "# Infraspecific name:", "")
                                strain = strline.replace(
                                    "strain=", "").strip()
                            elif gline.startswith('# Isolate:'):
                                isoline = gline.replace(
                                    "# Isolate:", "")
                                isolate = isoline.replace(
                                    "strain=", "").strip()
                            if isolate != '' and strain != '':
                                break

                wgsmaster_file = os.path.join(
                    line_split[1], genome_dir_id + '_wgsmaster.gbff')
                if os.path.exists(wgsmaster_file):
                    with open(wgsmaster_file, 'r') as wfile:
                        for wline in wfile:
                            if pattern_strain.match(wline):
                                wstrain = wline.replace(
                                    "/strain=", '').replace('"', '').rstrip().lstrip()
                            elif pattern_isolate.match(wline):
                                wsisolate = wline.replace(
                                    "/isolate=", '').replace('"', '').rstrip().lstrip()
                            if wsisolate != '' and wstrain != '':
                                break

                genomic_file = os.path.join(
                    line_split[1], genome_dir_id + '_genomic.gbff')
                if os.path.exists(genomic_file):
                    with open(genomic_file, 'r') as gfile:
                        for gline in gfile:
                            if pattern_strain.match(gline):
                                gstrain = gline.replace(
                                    "/strain=", '').replace('"', '').rstrip().lstrip()
                                gstrains.append(gstrain)
                            elif pattern_isolate.match(gline):
                                gisolate = gline.replace(
                                    "/isolate=", '').replace('"', '').rstrip().lstrip()
                                gisolates.append(gisolate)

                if 'substr.' in strain:
                    print "{1} strain {0}".format(strain, genome_id)
                    strain = strain.split("substr.")[0].rstrip()
                if 'substr.' in gstrain:
                    print "{1} gstrain {0}".format(gstrain, genome_id)
                    gstrain = gstrain.split("substr.")[0].rstrip()
                if 'substr.' in wstrain:
                    print "{1} wstrain {0}".format(wstrain, genome_id)
                    wstrain = wstrain.split("substr.")[0].rstrip()

                combined_strain = []
                #==============================================================
                # print strain
                # print isolate
                # print wstrain
                # print wsisolate
                # print gstrains
                # print gisolates
                #==============================================================

                combined_strain.extend(self.standardise_strain([strain]))
                combined_strain.extend(self.standardise_strain([isolate]))
                combined_strain.extend(self.standardise_strain([wstrain]))
                combined_strain.extend(self.standardise_strain([wsisolate]))
                combined_strain.extend(self.standardise_strain(gstrains))
                combined_strain.extend(self.standardise_strain(gisolates))

                stripped_combined_strain = map(str.strip, combined_strain)

                #==============================================================
                # print (genome_id, ';'.join(
                #     set(filter(None, stripped_combined_strain))))
                #==============================================================

                if genome_id.startswith('GCA'):
                    typemat = self.genbank_dictionary.get(genome_id)
                else:
                    typemat = self.refseq_dictionary.get(genome_id)

                outf.write("{0}\t{1}\t{2}\t{3}\n".format(
                    genome_id, species, ';'.join(set(filter(None, stripped_combined_strain))), typemat))
        outf.close()

    def standardise_strain(self, list_strain):
        results = []
        for rawst in list_strain:
            # clean strain
            if 'type strain:' in rawst:
                rawst = rawst.replace('type strain:', '')
            if '(=' in rawst and ')' in rawst:
                rawst = rawst.replace('(=', '=').replace(')', '=')

            # split strain using special character
            temp_results = []
            if '/' in rawst:
                temp_results.append(rawst)
                temp_results.extend(rawst.split("/"))
            else:
                temp_results = [rawst]
            for temp_rawst in temp_results:
                results.extend(re.split(';|,|=', temp_rawst))
        return results


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--genome_dir_file',
                        help='file indicating path to all genomes directories')
    parser.add_argument('--assembly_summary_genbank',
                        help='assembly_summary_genbank.txt downloaded from NCBI')
    parser.add_argument('--assembly_summary_refseq',
                        help='assembly_summary_refseq.txt downloaded from NCBI')

    parser.add_argument('--out_dir', help='output directory')

    args = parser.parse_args()

    try:
        p = StrainParser(args.assembly_summary_genbank,
                         args.assembly_summary_refseq)
        p.run(args.genome_dir_file, args.out_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
