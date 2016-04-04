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

__prog_name__ = 'update_database_from_ftp.py'
__prog_desc__ = ('Update the GTDB with the latest genome downloaded from FTP.' +
                 'Before this update, make sure all metadata have been generated and CheckM did run on all new genomes')

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
import argparse
import datetime

from dateutil.parser import parse
from database_configuration import GenomeDatabaseConnectionFTPUpdate


class UpdateGTDBDatabase(object):

    def __init__(self, db, date):

        self.db = db
        self.id_database = 3  # By default we set the id to genbank (it is either 2 or 3
        if db == "refseq":
            self.id_database = 2
        self.domains = ["archaea", "bacteria"]
        self.report_database_update = open(
            "report_beta_{0}_{1}_update_db.log".format(db, date), "w")

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate()
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

    def runUpdate(self, checkm, genome_dirs_file, dl_date):

        update_date = self.parse_date(dl_date)
        dict_existing_records = self._populateExistingRecords()
        list_checkm_records = self._populateNewRecords(checkm)
        genome_dirs_dict = self._populateGenomeDirs(genome_dirs_file)
        # Check if the genome is an existing genome
        # short_checkm_records list the records that are either to add or version
        short_checkm_records = self._updateExistingGenomes(dict_existing_records, list_checkm_records, genome_dirs_dict)
        self.temp_con.commit()
        self._addOrVersionNewGenomes(dict_existing_records, short_checkm_records, genome_dirs_dict, update_date)
        self.temp_con.commit()

        # Because we have added and updated script we repopulate dict_existing_records
        dict_existing_records = self._populateExistingRecords()
        self._checkPathorRemoveRecord(dict_existing_records, genome_dirs_dict, short_checkm_records)
        self.temp_con.commit()
        self.report_database_update.close()

    def _checkPathorRemoveRecord(self, dict_existing_records, genome_dirs_dict, list_checkm_records):
        count = 1
        for record in dict_existing_records:
            print "_checkPathorRemoveRecord: {0}/{1}".format(count, len(dict_existing_records))
            count += 1
            # if the record was part of the checkm file, it has already been updated
            if record in list_checkm_records:
                continue
            else:
                if record not in genome_dirs_dict:
                    self._removeRecord(record)
                else:
                    self._checkPathRecord(record, dict_existing_records[record], genome_dirs_dict[record])

    def _removeRecord(self, record):
        self.report_database_update.write(
            "{0}\t{1}\tremoved\n".format(self.db, record))
        query = ("SELECT gl.id,gl.name,us.username " +
                 "FROM genome_lists as gl " +
                 "LEFT JOIN genome_list_contents as glc on glc.list_id=gl.id " +
                 "LEFT JOIN users as us on us.id=gl.owner_id " +
                 "LEFT JOIN genomes as ge on ge.id=glc.genome_id " +
                 "WHERE ge.name like '{0}'".format(record))
        self.temp_cur.execute(query)
        raw_results = self.temp_cur.fetchall()
        self.report_database_update.write("{0}\n".format(record))
        if len(raw_results) > 0:
            self.report_database_update.write("modified list(s) [{0}]:\n".format(record))
            for result in raw_results:
                self.report_database_update.write(
                    "list_id:{0}\tlist_name:{1}\t,list_owner:{2}\n".format(*result))
            self.report_database_update.write("###########\n")
        else:
            self.report_database_update.write("No list has been modified\n")
            self.report_database_update.write("###########\n")
        query_delete = (
            "DELETE FROM genomes WHERE name LIKE '{0}'".format(record))
        self.temp_cur.execute(query_delete)

    def _checkPathRecord(self, record, path_in_db, path_in_folder):
        if path_in_db not in path_in_folder:
            path_in_folder = re.sub(r"(^.+\/)(archaea\/|bacteria\/)", r"\g<2>", path_in_folder)
            path_in_folder += "/" + os.path.basename(path_in_folder)
            path_in_db = re.sub(r"(.+)(_genomic.fna)", r"\g<1>", path_in_db)
            query = "update genomes set fasta_file_location = replace(fasta_file_location, '{0}', '{1}') where id_at_source like '{2}'".format(
                    path_in_db, path_in_folder, record)
            self.report_database_update.write("{0}\t{1}\tupdate path\t{2}\t{3}\n".format(self.db, record, path_in_db, path_in_folder))
            self.temp_cur.execute(query)
            query = "update genomes set genes_file_location = replace(genes_file_location, '{0}', '{1}') where id_at_source like '{2}'".format(
                    path_in_db, path_in_folder, record)
            self.temp_cur.execute(query)

    def _addOrVersionNewGenomes(self, dict_existing_records, list_checkm_records, genome_dirs_dict, update_date):
        count = 1
        for checkm_record in list_checkm_records:
            print "_addOrVersionNewGenomes: {0}/{1}".format(count, len(list_checkm_records))
            count += 1
            if (checkm_record not in dict_existing_records) and (checkm_record in genome_dirs_dict):
                check_record_base = checkm_record.rsplit(".", 1)[0]
                id_record = self._checkPreviousVersion(check_record_base)
                if id_record < 0:  # -1
                    # we add the genome to the database
                    self._addNewGenomes(checkm_record, genome_dirs_dict, update_date)
                else:
                    self._addNewGenomes(checkm_record, genome_dirs_dict, update_date, id_record)

    def _addNewGenomes(self, checkm_record, genome_dirs_dict, update_date, id_record=None):
        list_genome_details = [checkm_record]
        list_genome_details.append('')  # description
        list_genome_details.append(True)
        list_genome_details.append(None)
        fasta_file_path = os.path.join(genome_dirs_dict[checkm_record],
                                       os.path.basename(genome_dirs_dict[checkm_record]) + "_genomic.fna")
        fasta_file_path_shorten = re.sub(r"(.+/)(archaea\/|bacteria\/)", r"\g<2>", fasta_file_path)
        list_genome_details.append(fasta_file_path_shorten)
        list_genome_details.append(self.sha256Calculator(fasta_file_path))
        list_genome_details.append(self.id_database)
        list_genome_details.append(checkm_record)
        list_genome_details.append(update_date)
        list_genome_details.append(True)
        list_genome_details.append(update_date)
        gene_file_path = os.path.join(genome_dirs_dict[checkm_record],
                                      os.path.basename(genome_dirs_dict[checkm_record]) + "_protein.faa")
        gene_file_path_shorten = re.sub(r"(.+/)(archaea\/|bacteria\/)", r"\g<2>", gene_file_path)
        list_genome_details.append(gene_file_path_shorten)
        list_genome_details.append(self.sha256Calculator(gene_file_path))
        if id_record is None:
            self.report_database_update.write("{0}\t{1}\tadd\n".format(self.db, checkm_record))
            self.temp_cur.execute("INSERT INTO genomes " +
                                  "(name,description,owned_by_root,owner_id,fasta_file_location,fasta_file_sha256,genome_source_id,id_at_source,date_added,has_changed,last_update,genes_file_location,genes_file_sha256) " +
                                  "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", list_genome_details)
        else:
            self.report_database_update.write("{0}\t{1}\tversion\n".format(self.db, checkm_record))
            self.temp_cur.execute("UPDATE genomes " +
                                  "SET name = %s,description = %s, " +
                                  "owned_by_root = %s,owner_id = %s,  " +
                                  "fasta_file_location = %s,fasta_file_sha256 = %s, " +
                                  "genome_source_id = %s,id_at_source = %s,  " +
                                  "date_added = %s,has_changed = %s,  " +
                                  "last_update = %s,genes_file_location = %s,  " +
                                  "genes_file_sha256 = %s WHERE id = {0}; ".format(id_record), list_genome_details)

    def _updateExistingGenomes(self, dict_existing_records, list_checkm_records, genome_dirs_dict):
        new_list_checkm_records = []
        count = 1
        for checkm_record in list_checkm_records:
            print "_updateExistingGenomes: {0}/{1}".format(count, len(list_checkm_records))
            count += 1
            if (checkm_record in dict_existing_records) and (checkm_record in genome_dirs_dict):
                self.report_database_update.write("{0}\t{1}\tupdate protein file\n".format(self.db, checkm_record))
                path_gtdb = re.sub(r"(^.+\/)(archaea\/|bacteria\/)", r"\g<2>", genome_dirs_dict[checkm_record])
                path_gtdb += "/" + os.path.basename(path_gtdb)
                path_database = re.sub(r"(.+)(_genomic.fna)", r"\g<1>", dict_existing_records[checkm_record])
#                 if checkm_record == 'GCA_001242845.1':
#                     print path_gtdb
#                     print path_database
#                     print "update genomes set fasta_file_location = replace(fasta_file_location, '{0}', '{1}') where name like '{2}'".format(
#                         path_database, path_gtdb, checkm_record)
                # If the record is in a different folder , we need to change it's path in the database
                if path_database not in path_gtdb:
                    query = "update genomes set fasta_file_location = replace(fasta_file_location, '{0}', '{1}') where name like '{2}'".format(
                            path_database, path_gtdb, checkm_record)
                    self.temp_cur.execute(query)
                    query = "update genomes set genes_file_location = replace(genes_file_location, '{0}', '{1}') where name like '{2}'".format(
                            path_database, path_gtdb, checkm_record)
                    self.temp_cur.execute(query)

                # if the records is in the Checkm folder that means genomics and protein files have changed. We need to re write their sha256 values
                genomic_files = glob.glob(genome_dirs_dict[checkm_record] + "/*_genomic.fna")
                if len(genomic_files) == 1:
                    genomic_file = genomic_files[0]
                    new_md5_genomic = self.sha256Calculator(genomic_file)
                    query = "update genomes set fasta_file_sha256 = '{0}' where name like '{1}'".format(
                        new_md5_genomic, checkm_record)
                    self.temp_cur.execute(query)
                gene_files = glob.glob(genome_dirs_dict[checkm_record] + "/*_protein.faa")
                if len(gene_files) == 1:
                    gene_file = gene_files[0]
                    new_md5_gene = self.sha256Calculator(gene_file)
                    query = "update genomes set genes_file_sha256 = '{0}' where name like '{1}'".format(
                        new_md5_gene, checkm_record)
                    self.temp_cur.execute(query)
            else:
                new_list_checkm_records.append(checkm_record)
        return new_list_checkm_records

    def _populateGenomeDirs(self, genome_dirs_file):
        with open(genome_dirs_file, 'r') as gen_file:
            gen_dict = {gen_line.split("\t")[0]: gen_line.split("\t")[1].strip()
                        for gen_line in gen_file}
        return gen_dict

    def _populateNewRecords(self, checkm):
        list_result = []
        checkm_fh = open(checkm, "rb")
        checkm_fh.readline()
        for line in checkm_fh:
            full_name = line.split("\t")[0]
            name = full_name.split("_")[0] + "_" + full_name.split("_")[1]
            list_result.append(name)
        return list_result

    def _populateExistingRecords(self):
        self.temp_cur.execute("SELECT  gen.name,gen.fasta_file_location " +
                              "FROM genomes as gen " +
                              "WHERE gen.genome_source_id = {0} ;".format(self.id_database))
        dict_records = {key: value for (key, value) in self.temp_cur}
        return dict_records


# Tools

    def _checkPreviousVersion(self, checkm_record):
        self.temp_cur.execute("SELECT  gen.id " +
                              "FROM genomes as gen " +
                              "WHERE gen.id_at_source like '{0}.%' ;".format(checkm_record))
        list_result = [record for (record,) in self.temp_cur]
        if len(list_result) > 0:
            return list_result[0]
        else:
            return -1

    def sha256Calculator(self, file_path):

        try:
            filereader = open(file_path, "rb")
        except:
            raise Exception("Cannot open Fasta file: " + file_path)
        m = hashlib.sha256()
        for line in filereader:
            m.update(line)
        sha256_checksum = m.hexdigest()
        filereader.close()
        return sha256_checksum

    def parse_date(self, date_text):
        try:
            datetime.datetime.strptime(date_text, '%Y-%m-%d')
        except ValueError:
            raise ValueError("Incorrect data format, should be YYYY-MM-DD")
        return parse(date_text)

if __name__ == "__main__":
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--checkm_profile_new_genomes', dest="checkm",
                        required=True, help='path_to_checkm file')
    parser.add_argument('--genome_dirs_file', dest="genome_dirs_file",
                        required=True, help='genome_dirs file listing all REcords')
    parser.add_argument('--database', dest="db", required=True, choices=["refseq", "genbank"],
                        help='RefSeq or Genbank')
    parser.add_argument('--ftp_download_date', dest="date", required=True,
                        help='Date when the FTP download has been run (format YYYY-MM-DD)')

    args = parser.parse_args()

    try:
        update_mngr = UpdateGTDBDatabase(args.db, args.date)
        update_mngr.runUpdate(args.checkm, args.genome_dirs_file, args.date)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
