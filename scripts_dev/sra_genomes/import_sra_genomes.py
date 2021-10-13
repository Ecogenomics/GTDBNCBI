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

__prog_name__ = 'import_sra_genomes.py'
__prog_desc__ = ('Add externally processed genomes to GTDB ')

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2016'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@qfab.org'
__status__ = 'Development'

import os
import shutil
import hashlib
import re
import glob
import gzip
import sys
import datetime
import argparse


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

__prog_name__ = 'update_refseq_from_ftp_files.py'
__prog_desc__ = ('Update the Refseq folder with the latest genome downloaded from FTP.' +
                 'Before this update, make sure the genome folder is in RW mode to be able to delete or add genomes.' +
                 'We assume that all paths are in the form bacteria_name/latest_assembly_Version/GCA or GCF number')

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2016'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@qfab.org'
__status__ = 'Development'

import os
import shutil
import hashlib
import re
import glob
import gzip
import sys
import datetime
import argparse
import tempfile

from database_configuration import GenomeDatabaseConnectionSRAUpdate
from biolib.checksum import sha256


class UpdateSRA(object):

    def __init__(self):
        self.temp_con = GenomeDatabaseConnectionSRAUpdate.GenomeDatabaseConnectionSRAUpdate()
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

    def copytree(self, src, dst, symlinks=False, ignore=None):
        for item in os.listdir(src):
            s = os.path.join(src, item)
            d = os.path.join(dst, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, symlinks, ignore)
            else:
                shutil.copy2(s, d)

    def parsecheckm(self, sra_path):
        checkmfile = os.path.join(sra_path, "checkm_profile.tsv")
        checkm_dict = {}
        with open(checkmfile, 'r') as chfile:
            chfile.readline()
            for line in chfile:
                line_tab = line.split('\t')
                print(line_tab)
                checkm_dict[line_tab[0]] = {"checkm_completeness": line_tab[5], "checkm_contamination": line_tab[6], "checkm_marker_lineage": line_tab[1], "checkm_genome_count": line_tab[2], "checkm_marker_set_count": line_tab[4], "checkm_marker_count": line_tab[3], "checkm_strain_heterogeneity": line_tab[7]}
        return checkm_dict

    def run(self, sra_path):

        query = ("SELECT last_auto_id from genome_sources where id =1")
        self.temp_cur.execute(query)
        last_id = int(self.temp_cur.fetchone()[0])

        checkm_dict_original = self.parsecheckm(sra_path)

        print(os.getlogin())
        sra_dirs = os.listdir(sra_path)
        for sra_dir in sra_dirs:
            temp_path = tempfile.mkdtemp()
            sra_dir = os.path.join(sra_path, sra_dir)
            print(sra_dir)
            if os.path.isdir(sra_dir):
                bins_dir = os.path.join(sra_dir, 'bins')
                genomes_bin = os.listdir(bins_dir)
                for genome in genomes_bin:
                    if genome.endswith(".fna"):
                        genome_prefix = genome[:-4]
                        last_id += 1
                        temp_user_dir = os.path.join(temp_path, "U_" + str(last_id))
                        os.mkdir(temp_user_dir)
                        print(temp_user_dir)
                        print(genome_prefix)

                        shutil.copyfile(os.path.join(bins_dir, genome), os.path.join(temp_user_dir, "U_" + str(last_id) + "_genomic.fna"))

                        metadata_dir = os.path.join(sra_dir, 'metadata', genome_prefix)
                        self.copytree(metadata_dir, temp_user_dir)

                        prodigal_dir = os.path.join(temp_user_dir, 'prodigal')
                        os.mkdir(prodigal_dir)

                        for old_name in glob.glob(sra_dir + "/prodigal/" + genome_prefix + "_*"):
                            new_name = old_name.replace(sra_dir + "/prodigal", prodigal_dir)
                            new_name = new_name.replace(genome_prefix, "U_" + str(last_id))
                            shutil.copy(old_name, new_name)

                        for old_name in glob.glob(sra_dir + "/pfam/" + genome_prefix + "_*"):
                            new_name = old_name.replace(sra_dir + "/pfam", prodigal_dir)
                            new_name = new_name.replace(genome_prefix, "U_" + str(last_id))
                            shutil.copy(old_name, new_name)

                        for old_name in glob.glob(sra_dir + "/tigrfam/" + genome_prefix + "_*"):
                            new_name = old_name.replace(sra_dir + "/tigrfam", prodigal_dir)
                            new_name = new_name.replace(genome_prefix, "U_" + str(last_id))
                            shutil.copy(old_name, new_name)

                        list_genome_details = [genome_prefix]
                        list_genome_details.append(genome_prefix + " (21/07/2016)")  # description
                        list_genome_details.append(False)
                        list_genome_details.append(30)
                        fasta_file_path = os.path.join(os.path.join(os.getlogin(), genome_prefix, genome_prefix + "_genomic.fna"))
                        list_genome_details.append(fasta_file_path)
                        list_genome_details.append(sha256(os.path.join(temp_user_dir, "U_" + str(last_id) + "_genomic.fna")))
                        list_genome_details.append(1)
                        list_genome_details.append(str(last_id))
                        list_genome_details.append("21-07-2016")
                        list_genome_details.append(False)
                        list_genome_details.append("21-07-2016")
                        gene_file_path = os.path.join(os.getlogin(), genome_prefix, "prodigal", genome_prefix + "_protein.faa")
                        list_genome_details.append(gene_file_path)
                        list_genome_details.append(sha256(os.path.join(temp_user_dir, "prodigal", "U_" + str(last_id) + "_protein.faa")))

                        self.temp_cur.execute("INSERT INTO genomes " +
                                              "(name,description,owned_by_root,owner_id,fasta_file_location,fasta_file_sha256,genome_source_id,id_at_source,date_added,has_changed,last_update,genes_file_location,genes_file_sha256) " +
                                              "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s) RETURNING id ", list_genome_details)
                        (new_gid,) = self.temp_cur.fetchone()
                        self.temp_cur.execute(
                            "INSERT INTO metadata_nucleotide (id) VALUES ({0})".format(new_gid))
                        self.temp_cur.execute(
                            "INSERT INTO metadata_genes (id) VALUES ({0})".format(new_gid))
                        self.temp_cur.execute(
                            "INSERT INTO metadata_taxonomy (id) VALUES ({0})".format(new_gid))
                        self.temp_cur.execute(
                            "INSERT INTO metadata_ssu (id) VALUES ({0})".format(new_gid))

                        # insertion of metadata
                        with open(os.path.join(temp_user_dir, "metadata.genome_nt.tsv")) as metntf:
                            for line in metntf:
                                line_tab = line.strip().split()
                                self.temp_cur.execute("UPDATE metadata_nucleotide set {0}=%s WHERE id =%s ".format(line_tab[0]), (line_tab[1], new_gid))

                        # insertion of metadata
                        with open(os.path.join(temp_user_dir, "metadata.genome_gene.tsv")) as metgenef:
                            for line in metgenef:
                                line_tab = line.strip().split()
                                self.temp_cur.execute("UPDATE metadata_genes set {0}=%s WHERE id =%s ".format(line_tab[0]), (line_tab[1], new_gid))

                        for key, value in list(checkm_dict_original.get(genome_prefix).items()):
                            self.temp_cur.execute("UPDATE metadata_genes set {0}=%s WHERE id =%s ".format(key), (value, new_gid))

                        shutil.copytree(os.path.join(temp_user_dir), os.path.join("/srv/db/gtdb/genomes/user", os.getlogin(), "U_" + str(last_id)))

                        self.temp_con.commit()
            shutil.rmtree(temp_path)

if __name__ == "__main__":
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sra_folder', dest="sra_path", required=True,
                        help='path to the folder containing all SRA information')
    args = parser.parse_args()

    try:
        update_mngr = UpdateSRA()
        update_mngr.run(args.sra_path)

    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
