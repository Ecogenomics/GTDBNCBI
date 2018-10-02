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

__prog_name__ = 'select_gtdb_type_genomes.py'
__prog_desc__ = 'SELECT GTDB genomes as species reference for GTDB Tree'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2018'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

import os
import sys
import argparse
import tempfile
import re
import shutil

from biolib.taxonomy import Taxonomy


class TypeSelector(object):
    """Select GTDB type genomes as flag holder for species naming
      """

    def __init__(self):
        """Initialization."""
        pass
        #path to genomes
        self.path_to_genomes = {}
        with open("genome_dirs.tsv") as gdr:
            for line in gdr:
                infos = line.strip().split("\t")
                self.path_to_genomes[infos[0]] = os.path.join(infos[1], os.path.basename(infos[1]) + "_genomic.fna")

    def trunc_at(self, s, d, n=2):
        "Returns s truncated at the n'th (3rd by default) occurrence of the delimiter, d."
        return d.join(s.split(d, n)[:n])

    def parse_fastani_results(self, fastout_file):
        """ Parse the fastani output file

        Parameters
        ----------
        fastout_file : fastani output file.

        Returns
        -------
        dictionary
            dict_results[user_g]={"ref_genome":ref_genome,"ani":ani}
        """
        dict_results = {}
        with open(fastout_file) as fastfile:
            for line in fastfile:
                info = line.strip().split(" ")
                ref_genome = os.path.basename(info[1]).replace("_genomic.fna", "")
                ref_genome = self.trunc_at(ref_genome, '_')
                user_g = os.path.basename(info[0]).replace("_genomic.fna", "")
                user_g = self.trunc_at(user_g, '_')
                ani = float(info[2])
                if user_g in dict_results:
                    dict_results[user_g].update({ref_genome: ani})
                else:
                    dict_results[user_g] = {ref_genome: ani}

        return dict_results

    def calculate_fastani_distance(self, genome, temp_genome):
        """ Calculate the FastANI distance between all user genomes and the reference to classfy them at the species level

        Parameters
        ----------
        list_leaf : List of leaves uncluding one or many user genomes and one reference genome.
        genomes : Dictionary of user genomes d[genome_id] -> FASTA file

        Returns
        -------
        dictionary
            dict_results[user_g]={"ref_genome":ref_genome,"mash_dist":mash_dist}

        """
        #try:
        tmp_output_dir = tempfile.mkdtemp()
        query_list_file = open(os.path.join(tmp_output_dir, 'query_list.txt'),'w')
        ref_list_file = open(os.path.join(tmp_output_dir, 'ref_list.txt'),'w')

        query_list_file.write(('{0}\n'.format(self.path_to_genomes.get(genome[3:]))))
        ref_list_file.write(('{0}\n'.format(self.path_to_genomes.get(temp_genome[3:]))))

        ref_list_file.close()
        query_list_file.close()

        if not os.path.isfile(os.path.join(tmp_output_dir, 'query_list.txt')) or not os.path.isfile(os.path.join(tmp_output_dir, 'ref_list.txt')):
            raise

        cmd = 'fastANI --ql {0} --rl {1} -o {2} > /dev/null 2>{3}'.format(os.path.join(tmp_output_dir, 'query_list.txt'),
                                                                          os.path.join(tmp_output_dir, 'ref_list.txt'),
                                                                          os.path.join(tmp_output_dir, 'results.tab'),
                                                                          os.path.join(tmp_output_dir, 'error.log'))
        os.system(cmd)
        #print cmd

        if not os.path.isfile(os.path.join(tmp_output_dir, 'results.tab')):
            errstr = 'FastANI has stopped:\n'
            if os.path.isfile(os.path.join(tmp_output_dir, 'error.log')):
                with open(os.path.join(tmp_output_dir, 'error.log')) as debug:
                    for line in debug:
                        finalline = line
                    errstr += finalline
            raise ValueError(errstr)

        dict_parser_distance = self.parse_fastani_results(os.path.join(tmp_output_dir, 'results.tab'))
        #shutil.rmtree(tmp_output_dir)
        return dict_parser_distance

        #=======================================================================
        # except ValueError as error:
        #     if os.path.exists(tmp_output_dir):
        #         shutil.rmtree(tmp_output_dir)
        #     raise
        # except Exception as error:
        #     if os.path.exists(tmp_output_dir):
        #         shutil.rmtree(tmp_output_dir)
        #     raise
        #=======================================================================

    def run(self, ncbi_taxonomy, metadata_file):
        """Main function."""
        
        # output table summary
        summary_file = open("summary_table.tsv",'w')
        summary_file.write("species_name\tis_binomial\t#_type_strains\t#_subsp_type_strains\t")
        summary_file.write("#_ncbi_ref_genomes\t#_ncbi_rep_genomes\tresolving_case\ttype_strain_gids\ttype_subsp_strain_gids\tncbi_reference_genome\tncbi_representative_genome\n")

        # read taxonomy file
        taxonomy = Taxonomy().read(ncbi_taxonomy)
        spe_dict = {}
        for genome_id, taxstring in taxonomy.iteritems():
            if not genome_id.startswith("U_"):
                if taxstring[6] in spe_dict:
                    spe_dict.get(taxstring[6]).append(genome_id)
                else:
                    spe_dict[taxstring[6]] = [genome_id]

        type_dict = {}
        with open(metadata_file) as mf:
            header = mf.readline().strip().split('\t')
            gtdb_type_material_index = header.index('gtdb_type_material')
            gtdb_type_material_sources_index = header.index('gtdb_type_material_sources')
            straininfo_strains_index = header.index('straininfo_strains')
            lpsn_strains_index = header.index('lpsn_strains')
            dsmz_strains_index = header.index('dsmz_strains')
            ncbi_type_material_index = header.index('ncbi_type_material')
            ncbi_organism_name_index = header.index('ncbi_organism_name')
            ncbi_strain_identifiers_index = header.index('ncbi_strain_identifiers')
            ncbi_refseq_category_index = header.index('ncbi_refseq_category')
            checkm_completeness_index = header.index('checkm_completeness')
            checkm_contamination_index = header.index('checkm_contamination')


            for line in mf:
                line_split = line.strip().split('\t')
                gtdb_id = line_split[0]
                type_dict[gtdb_id] = {"gtm": line_split[gtdb_type_material_index],
                                      "gtms": line_split[gtdb_type_material_sources_index],
                                      "ntm": line_split[ncbi_type_material_index],
                                      "ss": line_split[straininfo_strains_index],
                                      'ls': line_split[lpsn_strains_index],
                                      'non': line_split[ncbi_organism_name_index],
                                      'nsi': line_split[ncbi_strain_identifiers_index],
                                      'ds': line_split[dsmz_strains_index],
                                      'nsc': line_split[ncbi_refseq_category_index],
                                      'comp': line_split[checkm_completeness_index],
                                      'cont': line_split[checkm_contamination_index]}

        onetype = 0
        multitype = 0
        zerotype = 0
        totalspe = 0
        zerotypemultisubsp = 0
        #aniflagfile = open("aniflagfile.txt",'w')
        binomialname = open("binnames.txt",'w')
        nonbinomialname = open("nobinnames.txt",'w')
        notypemultistrfile = open("notypemultistrfile.txt",'w')

        for spe, gids in spe_dict.iteritems():
            matchObj = re.search("^s__[a-zA-Z]+ [a-z]+$", spe)
            summary_list = [None]*11
            summary_list[0]=spe
            if matchObj:
                summary_list[1]="True"
                binomialname.write('{0}\n'.format(spe))
                totalspe += 1
                counttype = 0
                typespeciesgids = []
                typestrainsgids = []
                refgenomegids = []
                repgenomegids = []
                for gid in gids:
                    if (type_dict.get(gid).get('gtm') == 't' and
                        ' subsp. ' not in type_dict.get(gid).get('ls') and
                        ' subsp. ' not in type_dict.get(gid).get('ss') and
                        ' subsp. ' not in type_dict.get(gid).get('ds')):
                            counttype += 1
                            typespeciesgids.append(gid)
                    elif type_dict.get(gid).get('gtm') == 't':
                        typestrainsgids.append(gid)
                    elif type_dict.get(gid).get('nsc').startswith('reference'):
                        refgenomegids.append(gid)
                    elif type_dict.get(gid).get('nsc').startswith('representative'):
                        repgenomegids.append(gid)
                summary_list[2]=str(len(typespeciesgids))
                summary_list[3]=str(len(typestrainsgids))
                summary_list[4]=str(len(refgenomegids))
                summary_list[5]=str(len(repgenomegids))
                
                summary_list[7]=','.join(typespeciesgids)
                summary_list[8]=','.join(typestrainsgids)
                summary_list[9]=','.join(refgenomegids)
                summary_list[10]=','.join(repgenomegids)
                if counttype == 1:
                    onetype += 1
                    summary_list[6]='case A'
                elif counttype == 0:
                    zerocounttype = 0
                    for gid in gids:
                        if type_dict.get(gid).get('gtm') == 't':
                            zerocounttype += 1
                    if zerocounttype == 1:
                        onetype += 1
                        summary_list[6] = 'case C'
                    elif zerocounttype == 0:
                        zerotype += 1
                        #case E is processed later 
                        summary_list[6] = 'case E'
                    else:
                        multitype += 1
                        zerotypemultisubsp += 1
                        notypemultistrfile.write("{0}\n".format(spe))
                        summary_list[6] = 'case D'

                else:
                    aniissue = False
                    summaryanilist = []
                    conflictanilist = []

                    
                    for gid in typespeciesgids:
                        #=======================================================
                        # templist = list(typespeciesgids)
                        # templist.remove(gid)
                        # for tempgid in templist:
                        #     dict_gid = self.calculate_fastani_distance(gid, tempgid)
                        #     if gid[3:] in dict_gid and dict_gid.get(gid[3:]).get(tempgid[3:]) >= 95 :
                        #         summaryanilist.append("{0} {1}:{2}".format(gid, tempgid, dict_gid.get(gid[3:]).get(tempgid[3:])))
                        #     else:
                        #         aniissue = True
                        #         summaryanilist.append("{0} {1}: < 95%".format(gid, tempgid))
                        #         conflictanilist.append(gid)
                        #         conflictanilist.append(tempgid)
                        # if aniissue:
                        #     aniflagfile.write("################")
                        #     for item in summaryanilist:
                        #         aniflagfile.write('{0}\n'.format(item))
                        #     aniflagfile.write("################")
                        #=======================================================

                        if gid in ['GB_GCA_000615205.1','GB_GCA_000615205.1','GB_GCA_000720375.1','GB_GCA_000767215.1','GB_GCA_001313965.1','GB_GCA_001315225.1','GB_GCA_001315405.1','GB_GCA_001418385.1','GB_GCA_001940295.1','GB_GCA_001940295.1','GB_GCA_900156385.1','GB_GCA_900156385.1','GB_GCA_900227845.1','GB_GCA_900227845.1','RS_GCF_000014725.1','RS_GCF_000025985.1','RS_GCF_000025985.1','RS_GCF_000176795.1','RS_GCF_000187855.1','RS_GCF_000187855.1','RS_GCF_000197735.1','RS_GCF_000222765.1','RS_GCF_000222765.1','RS_GCF_000248355.1','RS_GCF_000258445.1','RS_GCF_000258445.1','RS_GCF_000260415.1','RS_GCF_000260415.1','RS_GCF_000260415.1','RS_GCF_000276905.1','RS_GCF_000276905.1','RS_GCF_000276905.1','RS_GCF_000300455.3','RS_GCF_000300455.3','RS_GCF_000369105.1','RS_GCF_000373025.1','RS_GCF_000382065.1','RS_GCF_000382065.1','RS_GCF_000382065.1','RS_GCF_000382065.1','RS_GCF_000382065.1','RS_GCF_000382065.1','RS_GCF_000383995.1','RS_GCF_000385925.1','RS_GCF_000385925.1','RS_GCF_000413475.1','RS_GCF_000413475.1','RS_GCF_000474095.1','RS_GCF_000474095.1','RS_GCF_000514775.1','RS_GCF_000597745.1','RS_GCF_000597745.1','RS_GCF_000597745.1','RS_GCF_000597745.1','RS_GCF_000597745.1','RS_GCF_000597745.1','RS_GCF_000612885.1','RS_GCF_000632805.1','RS_GCF_000719105.1','RS_GCF_000720555.1','RS_GCF_000721275.1','RS_GCF_000741335.1','RS_GCF_000741495.1','RS_GCF_000759445.1','RS_GCF_000771265.1','RS_GCF_000771405.1','RS_GCF_000803315.1','RS_GCF_000816845.1','RS_GCF_000829055.1','RS_GCF_000829055.1','RS_GCF_000950315.1','RS_GCF_001039055.1','RS_GCF_001302565.1','RS_GCF_001418405.1','RS_GCF_001431675.1','RS_GCF_001431675.1','RS_GCF_001431675.1','RS_GCF_001431675.1','RS_GCF_001431675.1','RS_GCF_001431675.1','RS_GCF_001433735.1','RS_GCF_001433735.1','RS_GCF_001436535.1','RS_GCF_001436535.1','RS_GCF_001436535.1','RS_GCF_001437155.1','RS_GCF_001437155.1','RS_GCF_001437155.1','RS_GCF_001507305.1','RS_GCF_001514065.1','RS_GCF_001550305.1','RS_GCF_001552375.1','RS_GCF_001558415.1','RS_GCF_001558415.1','RS_GCF_001570385.1','RS_GCF_001570625.1','RS_GCF_001591205.1','RS_GCF_001591205.1','RS_GCF_001591205.1','RS_GCF_001591205.1','RS_GCF_001591205.1','RS_GCF_001591205.1','RS_GCF_001591785.1','RS_GCF_001591785.1','RS_GCF_001598835.1','RS_GCF_001598835.1','RS_GCF_001636475.1','RS_GCF_001636475.1','RS_GCF_001636495.1','RS_GCF_001648575.1','RS_GCF_001729765.1','RS_GCF_001729765.1','RS_GCF_001890655.1','RS_GCF_001890655.1','RS_GCF_001894925.1','RS_GCF_001997185.1','RS_GCF_001997185.1','RS_GCF_001997185.1','RS_GCF_001997185.1','RS_GCF_001997185.1','RS_GCF_001997185.1','RS_GCF_001997325.1','RS_GCF_001997325.1','RS_GCF_002102225.1','RS_GCF_002104765.1','RS_GCF_002217195.1','RS_GCF_002217245.1','RS_GCF_002217275.1','RS_GCF_002243645.1','RS_GCF_002243645.1','RS_GCF_002335735.1','RS_GCF_002530675.1','RS_GCF_900078675.2','RS_GCF_900100425.1','RS_GCF_900100995.1','RS_GCF_900100995.1','RS_GCF_900107365.1','RS_GCF_900109485.1','RS_GCF_900109485.1','RS_GCF_900110015.1','RS_GCF_900110795.1','RS_GCF_900113425.1','RS_GCF_900113425.1','RS_GCF_900113425.1','RS_GCF_900113425.1','RS_GCF_900113425.1','RS_GCF_900113425.1','RS_GCF_900129595.1','RS_GCF_900186865.1','RS_GCF_900186865.1','RS_GCF_900186865.1','RS_GCF_900186865.1','RS_GCF_900186865.1','RS_GCF_900186865.1']:
                        #if gid in conflictanilist:
                            summary_list[6] = 'case B-1'
                            break

                    if len(refgenomegids) == 1 or len(repgenomegids) == 1 and summary_list[6] != 'case B-1':
                        onetype += 1
                        summary_list[6]='case B-2'
                    elif summary_list[6] != 'case B-1':
                        selected_gid=''
                        bestquality = -1000
                        for gid in typespeciesgids:
                            quality = float(type_dict.get(gid).get('comp'))-5*float(type_dict.get(gid).get('cont'))
                            if quality > bestquality:
                                bestquality = quality
                                selected_gid = gid
                        multitype += 1
                        summary_list[6] = 'case B-3'
                        
                if summary_list[6] == 'case E':
                    if len(refgenomegids) == 1 or len(repgenomegids) == 1 :
                        summary_list[6] = 'case E-1'
                    else:
                        summary_list[6] = 'case E-2'
                    

            else:
                summary_list[1] = "False"
                typespeciesgids = []
                typestrainsgids = []
                nonbinomialname.write('{0}\n'.format(spe))
                for gid in gids:
                    if (type_dict.get(gid).get('gtm') == 't' and
                        ' subsp. ' not in type_dict.get(gid).get('ls') and
                        ' subsp. ' not in type_dict.get(gid).get('ss') and
                        ' subsp. ' not in type_dict.get(gid).get('ds')):
                            counttype += 1
                            typespeciesgids.append(gid)
                    elif type_dict.get(gid).get('gtm') == 't':
                        typestrainsgids.append(gid)
                summary_list[2]=str(len(typespeciesgids))
                summary_list[3]=str(len(typestrainsgids))
                summary_list[4]=str(len(refgenomegids))
                summary_list[5]=str(len(repgenomegids))
                
                summary_list[7]=','.join(typespeciesgids)
                summary_list[8]=','.join(typestrainsgids)
                summary_list[9]=','.join(refgenomegids)
                summary_list[10]=','.join(repgenomegids)
            summary_file.write("{0}\n".format("\t".join(item or '' for item in summary_list)))
        #aniflagfile.close()
        summary_file.close()
        print "one type: {0}".format(onetype)
        print "zero type: {0}".format(zerotype)
        print "multiple types: {0}".format(multitype)
        print "total species: {0}".format(totalspe)
        print "\n\nzero type multi strains: {0}".format(zerotypemultisubsp)


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ncbi_taxonomy', help='TSV file of the NCBI Taxonomy exported from GTDB.')
    parser.add_argument('--metadata_file', help='Metadata file from GTDB')

    args = parser.parse_args()

    try:
        typeselector = TypeSelector()
        typeselector.run(args.ncbi_taxonomy, args.metadata_file)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
