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

__prog_name__ = 'update_gtdb_folder.py'
__prog_desc__ = ('Update the gtdb folder with the latest genome downloaded from FTP.' +
                 'Before this update, make sure the genome folder is in RW mode to be able to delete or add genomes.' +
                 'We assume that all paths are in the form bacteria_name/latest_assembly_Version/GCA or GCF number')

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2015'
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

from dateutil.parser import parse
from gtdblite import GenomeDatabaseConnection


class UpdateRefSeqFolder(object):

    def __init__(self):

        self.temp_con = GenomeDatabaseConnection.GenomeDatabaseConnection()
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

        self.domains = ["archaea", "bacteria"]

        self.genomic_ext = "_genomic.fna.gz"
        self.protein_ext = "_protein.faa.gz"
        self.fastaExts = (self.genomic_ext, self.protein_ext)
        self.extensions = ("_feature_table.txt.gz", "_genomic.gbff.gz",
                           "_genomic.gff.gz", "_protein.gpff.gz", "_wgsmaster.gbff.gz")
        self.reports = ("_assembly_report.txt", "_assembly_stats.txt")
        self.allExts = self.fastaExts + self.extensions + self.reports
        self.allbutFasta = self.extensions + self.reports

        self.report_gcf = open("report_gcf.log", "w")
        self.stats_update = open("stats_update.log", "w")

    def runComparison(self, ftp_refseq, new_refseq, ftp_genome_dirs, old_genome_dirs, download_date):
        '''
        runComparison function is walking across all directories recursively
        only folder containing latest_assembly_versions but not containing _assembly_structure
        are of interest
        '''
        update_date = self.parse_date(download_date)

        # for each domain
        for domain in self.domains:

            # old_dict lists all records from the previous gtdb update
            with open(old_genome_dirs, 'r') as old_file:
                old_dict = {old_line.split("\t")[0]: old_line.split("\t")[1].strip()
                            for old_line in old_file if "/{0}/".format(self.domain) in old_line.split("\t")[1]}

            # new list list all records from the ftp folder and considered as
            # latest
            ftp_assembly_summary = os.path.join(
                ftp_refseq, domain, "assembly_summary.txt")
            with open(ftp_assembly_summary, 'r') as ftp_assembly_summary_file:
                ftp_assembly_summary_file.readline()
                new_list = [new_line.split(
                    "\t")[0] for new_line in ftp_assembly_summary_file if new_line.split("\t")[10] == "latest"]
            # new dict lists all records from FTP which are in new_list
            with open(ftp_genome_dirs, 'r') as new_genome_dirs_file:
                new_dict = {new_line.split("\t")[0]: new_line.split("\t")[1].strip()
                            for new_line in new_genome_dirs_file if "/{0}/".format(self.domain) in new_line.split("\t")[1] and new_line.split("\t")[0] in new_list}

            # new genomes in FTP
            added_dict = {added_key: new_dict[added_key] for added_key in list(
                set(new_dict.keys()) - set(old_dict.keys()))}
            self.addGenomes(added_dict)

            # delete genomes from the Database
            removed_dict = {removed_key: old_dict[removed_key] for removed_key in list(
                set(old_dict.keys()) - set(new_dict.keys()))}
            self.removeGenomes(removed_dict)

            intersect_list = list(
                set(old_dict.keys()).intersection(set(new_dict.keys())))
            self.compareGenomes(intersect_list, old_dict, new_dict)

    def addGenomes(self, added_dict, ftp_refseq, new_refseq, update_date):

        for gcf_record in added_dict:
            target_dir = added_dict[gcf_record].replace(ftp_refseq, new_refseq).replace(
                "/latest_assembly_versions", "")
            shutil.copytree(added_dict[
                            gcf_record], target_dir, ignore=shutil.ignore_patterns("*_assembly_structure"))
            self.report_gcf.write(
                "{0}\t{1}\tnew\n".format(self.domain.upper(), gcf_record))
            for compressed_file in glob.glob(target_dir + "/*.gz"):
                if os.path.isdir(compressed_file) == False:
                    inF = gzip.open(compressed_file, 'rb')
                    outF = open(compressed_file.replace(".gz", ""), 'wb')
                    outF.write(inF.read())
                    inF.close()
                    outF.close()
                    os.remove(compressed_file)
            # write all info for genomes
            list_genome_details = [gcf_record]
            list_genome_details.append('')  # description
            list_genome_details.append(True)
            list_genome_details.append(None)
            fasta_file_path = os.path.join(
                target_dir, os.path.basename(target_dir) + "_genomic.fna")

            fasta_file_path_shorten = re.sub(
                "^.+\/{0}\/".format(self.domain), "{0}/".format(self.domain), fasta_file_path)
            list_genome_details.append(fasta_file_path_shorten)
            list_genome_details.append(self.sha256Calculator(fasta_file_path))
            list_genome_details.append(2)
            list_genome_details.append(gcf_record)
            list_genome_details.append(update_date)
            list_genome_details.append(True)
            list_genome_details.append(update_date)
            gene_file_path = os.path.join(
                target_dir, os.path.basename(target_dir) + "_protein.faa")
            gene_file_path_shorten = re.sub(
                "^.+\/{0}\/".format(self.domain), "{0}/".format(self.domain), gene_file_path)
            list_genome_details.append(gene_file_path_shorten)
            list_genome_details.append(self.sha256Calculator(gene_file_path))
            query = "INSERT INTO genomes (name,description,owned_by_root,owner_id,fasta_file_location,fasta_file_sha256,genome_source_id,id_at_source,date_added,has_changed,last_update,genes_file_location,genes_file_sha256) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
            self.temp_cur.execute(query, list_genome_details)
            self.temp_con.commit()

    def removeGenomes(self, removed_dict):
        for gcf_record in removed_dict:
            self.report_gcf.write(
                "{0}\t{1}\tremoved\n".format(self.domain.upper(), gcf_record))
            query = ("SELECT gl.id,gl.name,us.username " +
                     "FROM genome_lists as gl " +
                     "LEFT JOIN genome_list_contents as glc on glc.list_id=gl.id " +
                     "LEFT JOIN users as us on us.id=gl.owner_id " +
                     "LEFT JOIN genomes as ge on ge.id=glc.genome_id " +
                     "WHERE ge.name like '{0}'".format(gcf_record))
            self.temp_cur.execute(query)
            raw_results = self.temp_cur.fetchall()
            self.stats_update.write("{0}\n".format(gcf_record))
            if len(raw_results) > 0:
                self.stats_update.write("modified list(s):\n")
                for result in raw_results:
                    self.stats_update.write(
                        "list_id:{0}\tlist_name:{1}\t,list_owner:{2}\n".format(*result))
            else:
                self.stats_update.write("No list has been modified\n")
            query_delete = (
                "delete from genomes where name like '{0}'".format(gcf_record))
            self.temp_cur.execute(query_delete)
            self.temp_con.commit()

    def compareGenomes(self, intersect_list, old_dict, new_dict):
        for gcf_record in intersect_list:
            gtdb_dir = old_dict.get(gcf_record)
            ftp_dir = new_dict.get(gcf_record)
            target_dir = ftp_dir.replace(self.origin, self.target).replace(
                "latest_assembly_versions/", "")
            self.readmd5Checksum(gtdb_dir, ftp_dir, target_dir, gcf_record)

    def readmd5Checksum(self, gtdb_dir, ftp_dir, target_dir, gcf_record):
        '''
        Compare the checksum of the file listed in the checksums.txt
        '''
        pathftpmd5 = os.path.join(ftp_dir, "md5checksums.txt")
        pathgtdbmd5 = os.path.join(gtdb_dir, "md5checksums.txt")
        target_pathnewmd5 = os.path.join(target_dir, "md5checksums.txt")
        status = []

        path_gtdb = re.sub(
            "^.+\/{0}\/".format(self.domain), "{0}/".format(self.domain), gtdb_dir)
        path_ftp = re.sub(
            "^.+\/{0}\/".format(self.domain), "{0}/".format(self.domain), target_dir)
        if not path_ftp in path_gtdb:
            status.append("moved")
            query = "update genomes set fasta_file_location = replace(fasta_file_location, '{0}', '{1}') where name like '{2}'".format(
                path_gtdb, path_ftp, gcf_record)
            query = "update genomes set genes_file_location = replace(genes_file_location, '{0}', '{1}') where name like '{2}'".format(
                path_gtdb, path_ftp, gcf_record)
            self.temp_cur.execute(query)
            self.temp_con.commit()

        ftpdict, ftpdict_fasta = self.parse_checksum(pathftpmd5)
        gtdbdict, gtdbdict_fasta = self.parse_checksum(pathgtdbmd5)

        if len(list(set(ftpdict_fasta.keys()).symmetric_difference(set(gtdbdict_fasta.keys())))) > 0:
            status.append("incomplete")
            return False
        else:
            ftp_folder = True
            for key, value in ftpdict_fasta.iteritems():
                if value == gtdbdict_fasta.get(key):
                    ftp_folder = False
                print gcf_record + "\t" + key + "\t" + gtdbdict_fasta.get(key)
            if ftp_folder:
                shutil.copytree(
                    ftp_dir, target_dir,
                    ignore=shutil.ignore_patterns("*_assembly_structure"))
                status.append("modified")
                for compressed_file in glob.glob(target_dir + "/*.gz"):
                    if os.path.isdir(compressed_file) == False:
                        inF = gzip.open(compressed_file, 'rb')
                        outF = open(compressed_file.replace(".gz", ""), 'wb')
                        outF.write(inF.read())
                        inF.close()
                        outF.close()
                        # if the protein file or the seq file change we need to
                        #    recalculate the md5 of those new files
                        if self.genomic_ext in compressed_file:
                            new_md5 = self.sha256Calculator(
                                compressed_file.replace(".gz", ""))
                            query = "update genomes set fasta_file_sha256 = '{0}' where name like '{1}'".format(
                                new_md5, gcf_record)
                        if self.protein_ext in compressed_file:
                            new_md5 = self.sha256Calculator(
                                compressed_file.replace(".gz", ""))
                            query = "update genomes set genes_file_sha256 = '{0}' where name like '{1}'".format(
                                new_md5, gcf_record)
                        os.remove(compressed_file)
                self.temp_cur.execute(query)
                self.temp_con.commit()

            else:
                shutil.copytree(
                    gtdb_dir, target_dir,
                    ignore=shutil.ignore_patterns("*_assembly_structure"))
                status.append("unmodified")

                # print "{0} to {1}".format(gtdb_dir, target_dir)
                checksum_changed = False
                for key, value in ftpdict.iteritems():
                    if value != gtdbdict.get(key):
                        checksum_changed = True

                        shutil.copy2(
                            os.path.join(ftp_dir, key), os.path.join(target_dir, key))
                        if key.endswith(".gz"):
                            inF = gzip.open(os.path.join(ftp_dir, key), 'rb')
                            try:
                                print os.path.join(target_dir, key).replace(".gz", "")
                                outF = open(
                                    os.path.join(target_dir, key).replace(".gz", ""), 'wb')
                            except IOError:
                                os.chmod(
                                    os.path.join(target_dir, key).replace(".gz", ""), 0775)
                                outF = open(
                                    os.path.join(target_dir, key).replace(".gz", ""), 'wb')
                            outF.write(inF.read())
                            inF.close()
                            outF.close()
                            os.remove(os.path.join(target_dir, key))
                        status.append("new_metadata")
                if checksum_changed:
                    try:
                        shutil.copy2(pathgtdbmd5, target_pathnewmd5)
                    except IOError:
                        os.chmod(target_pathnewmd5, 0775)
                        shutil.copy2(pathftpmd5, target_pathnewmd5)
                for report in self.reports:
                    target_files = glob.glob(
                        os.path.join(gtdb_dir, "*" + report))
                    ftp_files = glob.glob(os.path.join(ftp_dir, "*" + report))
                    if len(target_files) == 1 and len(ftp_files) == 1:
                        status = self.comparesha256(
                            ftp_files[0], target_files[0], status)
                    else:
                        print target_files
                        print ftp_files
                        print "IT SHOULDN'T HAPPEN"
        self.report_gcf.write("{0}\t{1}\t{2}\n".format(
            self.domain.upper(), gcf_record, ';'.join([x for x in set(status)])))

# Tools

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

    def comparesha256(self, ftp_file, target_file, status):
        original_checksum = hashlib.md5(
            open(ftp_file, 'rb').read()).hexdigest()
        gtdb_checksum = hashlib.md5(open(target_file, 'rb').read()).hexdigest()
        if original_checksum != gtdb_checksum:
            try:
                shutil.copy2(ftp_file, target_file)
            except IOError:
                os.chmod(target_file, 0775)
                shutil.copy2(ftp_file, target_file)
            status.append("new_metadata")
        return status

    def parse_checksum(self, md5File):
        out_dict, out_dict_fasta = {}, {}

        with open(md5File) as f:
            for line in f:
                split_line = line.rstrip().split("  ")
                header = split_line[1].replace("./", "")
                chksum = split_line[0]
                if header.endswith(self.fastaExts):
                    out_dict_fasta[header] = chksum
                if header.endswith(self.allbutFasta):
                    out_dict[header] = chksum
        return (out_dict, out_dict_fasta)

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
    parser.add_argument('--ftp_refseq_directory', dest="ftp_refseq", required=True,
                        help='base directory leading the the FTP repository for refseq')
    parser.add_argument('--new_refseq_directory', dest="new_refseq",
                        required=True, help='base directory leading the new repository for refseq')
    parser.add_argument('--ftp_genome_dirs_file', dest="ftp_genome_dirs", required=True,
                        help='metadata file listing all directories for the FTP folder (generated by genome_dirs.py)')
    parser.add_argument('--old_genome_dirs_file', dest="old_genome_dirs", required=True,
                        help='metadata file listing all directories from the previous NCBI update date  (generated by genome_dirs.py)')
    parser.add_argument('--ftp_download_date', dest="download_date", required=True,
                        help='Date when the FTP download has been run (format YYYY-MM-DD')

    args = parser.parse_args()

    try:
        update_mngr = UpdateRefSeqFolder()
        update_mngr.runComparison(
            args.ftp_refseq, args.new_refseq, args.ftp_genome_dirs, args.old_genome_dirs, args.download_date)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
