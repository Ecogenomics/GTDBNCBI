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

__prog_name__ = 'update_refseq_from_ftp_files_uncompressed.py'
__prog_desc__ = ('Update the Refseq folder with the latest genome downloaded from FTP.' +
                 'This is an updated version of update_refseq_from_ftp_files where the checksum comparison is made on unarchived files.' +
                 'Before this update, make sure the genome folder is in RW mode to be able to delete or add genomes.' +
                 'We assume that all paths are in the form bacteria_name/latest_assembly_Version/GCA or GCF number.')

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


class UpdateRefSeqFolder(object):

    def __init__(self, new_refseq_folder):

        self.domains = ["archaea", "bacteria"]
        #self.domains = ["bacteria"]

        self.genomic_ext = "_genomic.fna"
        self.protein_ext = "_protein.faa"
        self.cds_ext = "_cds_from_genomic.fna"
        self.rna_ext = "_rna_from_genomic.fna"

        self.fastaExts = (self.genomic_ext, self.protein_ext)
        self.extrafastaExts = (self.cds_ext, self.rna_ext)

        self.extensions = ("_feature_table.txt", "_genomic.gbff",
                           "_genomic.gff", "_protein.gpff", "_wgsmaster.gbff")
        self.reports = ("_assembly_report.txt", "_assembly_stats.txt", "_hashes.txt")
        self.allExts = self.fastaExts + self.extensions + self.reports
        self.allbutFasta = self.extensions + self.reports

        self.report_gcf = open(os.path.join(new_refseq_folder, "report_gcf.log"), "w", 1)
        self.stats_update = open(os.path.join(new_refseq_folder, "stats_update.log"), "w")

    def runComparison(self, ftp_refseq, new_refseq, ftp_genome_dirs, old_genome_dirs):
        '''
        runComparison function is walking across all directories recursively
        only folder containing latest_assembly_versions but not containing _assembly_structure
        are of interest
        '''

        # for each domain
        for domain in self.domains:

            # old_dict lists all records from the previous gtdb update
            with open(old_genome_dirs, 'r') as old_file:
                old_dict = {old_line.split("\t")[0]: old_line.split("\t")[1].strip()
                            for old_line in old_file if "/{0}/".format(domain) in old_line.split("\t")[1]}

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
                            for new_line in new_genome_dirs_file if "/{0}/".format(domain) in new_line.split("\t")[1] and new_line.split("\t")[0] in new_list}

            # new genomes in FTP
            added_dict = {added_key: new_dict[added_key] for added_key in list(
                set(new_dict.keys()) - set(old_dict.keys()))}
            self.addGenomes(added_dict, ftp_refseq, new_refseq, domain)

            # delete genomes from the Database
            removed_dict = {removed_key: old_dict[removed_key] for removed_key in list(
                set(old_dict.keys()) - set(new_dict.keys()))}
            self.removeGenomes(removed_dict, domain)

            intersect_list = list(
                set(old_dict.keys()).intersection(set(new_dict.keys())))
            self.compareGenomes(
                intersect_list, old_dict, new_dict, ftp_refseq, new_refseq, domain)

    def addGenomes(self, added_dict, ftp_refseq, new_refseq, domain):
        '''
        addGenomes function insert new genomes in the GTDB database. New genomes are present in the FTP folder
        but not in the previous version of GTDB.

        :TODO: Check if the new genome is a new version of an existing genome. in that case we overwrite the previous one
        and keep the same database id
        This will cause a conflict with the removeGenomes function.


        :param added_dict: dictionary of genomes to be added (genome_id:path to genome)
        :param ftp_refseq: base directory leading the the FTP repository for refseq
        :param new_refseq:base directory leading the new repository for refseq
        :param update_date:Date when the FTP download has been run (format YYYY-MM-DD)
        '''

        for gcf_record in added_dict:
            target_dir = added_dict[gcf_record].replace(ftp_refseq, new_refseq).replace(
                "/latest_assembly_versions", "")
            shutil.copytree(added_dict[
                            gcf_record], target_dir, symlinks=True, ignore=shutil.ignore_patterns("*_assembly_structure"))
            self.report_gcf.write(
                "{0}\t{1}\tnew\n".format(domain.upper(), gcf_record))
            for compressed_file in glob.glob(target_dir + "/*.gz"):
                if os.path.isdir(compressed_file) == False:
                    inF = gzip.open(compressed_file, 'rb')
                    outF = open(
                        self.rreplace(compressed_file, ".gz", "", 1), 'wb')
                    outF.write(inF.read())
                    inF.close()
                    outF.close()
                    os.remove(compressed_file)

    def removeGenomes(self, removed_dict, domain):
        '''
        removeGenomes function removes all outdated genomes from the gtdb database
        In addition it tracks the lists(name and owner) that have been modified while deleting those genomes

        :param removed_dict: dictionary of genomes to delete
        '''

        for gcf_record in removed_dict:
            self.report_gcf.write(
                "{0}\t{1}\tremoved\n".format(domain.upper(), gcf_record))

    def compareGenomes(self, intersect_list, old_dict, new_dict, ftp_refseq, new_refseq, domain):
        '''
        compare the genomes existing in both folders ( FTP folder and previous gtdb update).

        :param intersect_list:
        :param old_dict:
        :param new_dict:
        '''
        for gcf_record in intersect_list:
            gtdb_dir = old_dict.get(gcf_record)
            ftp_dir = new_dict.get(gcf_record)
            target_dir = ftp_dir.replace(ftp_refseq, new_refseq).replace(
                "latest_assembly_versions/", "")
            self.readmd5Checksum(
                gtdb_dir, ftp_dir, target_dir, gcf_record, domain)

    def readmd5Checksum(self, gtdb_dir, ftp_dir, target_dir, gcf_record, domain):
        '''
        Compare the checksum of the file listed in the checksums.txt
        '''
        pathftpmd5 = os.path.join(ftp_dir, "md5checksums.txt")
        pathgtdbmd5 = os.path.join(gtdb_dir, "md5checksums.txt")
        target_pathnewmd5 = os.path.join(target_dir, "md5checksums.txt")
        status = []

        tmp_ftp_dir = tempfile.mkdtemp()
        tmp_target = os.path.join(tmp_ftp_dir, os.path.basename(target_dir))
        shutil.copytree(ftp_dir, tmp_target, symlinks=True, ignore=shutil.ignore_patterns("*_assembly_structure"))
        for compressed_file in glob.glob(tmp_target + "/*.gz"):
            if os.path.isdir(compressed_file) == False:
                inF = gzip.open(compressed_file, 'rb')
                try:
                    outF = open(
                        self.rreplace(compressed_file, ".gz", "", 1), 'wb')
                except IOError:
                    os.chmod(
                        self.rreplace(compressed_file, ".gz", "", 1), 0o775)
                    outF = open(
                        self.rreplace(compressed_file, ".gz", "", 1), 'wb')
                outF.write(inF.read())
                inF.close()
                outF.close()
                os.remove(compressed_file)

        ftpdict, ftpdict_fasta, ftpdict_faa, ftpdict_extra_fasta = self.parse_checksum(tmp_target)
        gtdbdict, gtdbdict_fasta, gtdbdict_faa, gtdbdict_extra_fasta = self.parse_checksum(gtdb_dir)

        # if the genomic.fna.gz or the protein.faa.gz are missing, we set this
        # record as incomplete
        if len(list(set(ftpdict_fasta.keys()).symmetric_difference(set(gtdbdict_fasta.keys())))) > 0:
            print tmp_target
            print ftpdict_fasta.keys()
            print gtdb_dir
            print gtdbdict_fasta.keys()
            status.append("incomplete")
            shutil.copytree(tmp_target, target_dir, symlinks=True, ignore=shutil.ignore_patterns("*_assembly_structure"))
            # we unzip of gz file

        else:
            ftp_folder = False
            # check if genomic.fna.gz and protein.faa.gz are similar between
            # previous gtdb and ftp
            for key, value in ftpdict_fasta.iteritems():
                if value != gtdbdict_fasta.get(key):
                    ftp_folder = True
            # if one of the 2 files is different than the previous version , we
            # use the ftp record over the previous gtdb one , we then need to
            # re run the metadata generation
            if ftp_folder:
                shutil.copytree(
                    tmp_target, target_dir, symlinks=True,
                    ignore=shutil.ignore_patterns("*_assembly_structure"))
                status.append("modified")

            else:
                # The 2 main fasta files haven't changed so we can copy the old
                # gtdb folder over
                shutil.copytree(
                    gtdb_dir, target_dir, symlinks=True,
                    ignore=shutil.ignore_patterns("*_assembly_structure"))
                status.append("unmodified")

                # We check if all other file of this folder are the same.
                checksum_changed = False

                for key, value in ftpdict_faa.iteritems():
                    if value != gtdbdict_faa.get(key):
                        checksum_changed = True
                        shutil.copy2(
                            os.path.join(tmp_target, key), os.path.join(target_dir, key))
                        status.append("new_protein")

                for key, value in ftpdict.iteritems():
                    if value != gtdbdict.get(key):
                        checksum_changed = True
                        shutil.copy2(
                            os.path.join(tmp_target, key), os.path.join(target_dir, key))
                        status.append("new_metadata")

                for key, value in ftpdict_extra_fasta.iteritems():
                    if value != gtdbdict_extra_fasta.get(key):
                        checksum_changed = True
                        shutil.copy2(
                            os.path.join(tmp_target, key), os.path.join(target_dir, key))
                        status.append("new_cds_rna")
                # we copy the new checksum
                if checksum_changed:
                    try:
                        shutil.copy2(pathgtdbmd5, target_pathnewmd5)
                    except IOError:
                        os.chmod(target_pathnewmd5, 0o664)
                        shutil.copy2(pathftpmd5, target_pathnewmd5)
                for report in self.reports:
                    target_files = glob.glob(
                        os.path.join(target_dir, "*" + report))
                    ftp_files = glob.glob(os.path.join(ftp_dir, "*" + report))
                    if len(target_files) == 1 and len(ftp_files) == 1:
                        status = self.comparesha256(
                            ftp_files[0], target_files[0], status)
                    else:
                        print "########"
                        print target_dir
                        print target_files
                        print ftp_dir
                        print ftp_files
                        print "IT SHOULDN'T HAPPEN"
                        print "########"
        self.report_gcf.write("{0}\t{1}\t{2}\n".format(
            domain.upper(), gcf_record, ';'.join([x for x in set(status)])))
        shutil.rmtree(tmp_ftp_dir)

# Tools

    def rreplace(self, s, old, new, occurrence):
        '''
         Instead of using the normal replace function, we need to implement our own.
         Some folder are named with a .gz in the middle so we only need to replace the last .gz in the string name
        :param s:
        :param old:
        :param new:
        :param occurrence:
        '''
        li = s.rsplit(old, occurrence)
        return new.join(li)

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
        '''
        comparesha256 compares the report file
        :param ftp_file:
        :param target_file:
        :param status:

        '''
        original_checksum = hashlib.md5(
            open(ftp_file, 'rb').read()).hexdigest()
        gtdb_checksum = hashlib.md5(open(target_file, 'rb').read()).hexdigest()
        if original_checksum != gtdb_checksum:
            try:
                shutil.copy2(ftp_file, target_file)
            except IOError:
                os.chmod(target_file, 0o664)
                shutil.copy2(ftp_file, target_file)
            status.append("new_metadata")
        return status

    def parse_checksum(self, pathtodir):
        '''
        parse_checksum function parses the md5 checksum file.
        It returns 2 dictionaries {file:size} : one for the fna and faa files, one for the genbank files
        :param md5File:
        '''

        out_dict, out_dict_fasta, out_dict_faa, out_dict_extra_fasta = {}, {}, {}, {}

        for name in glob.glob(os.path.join(pathtodir, '*')):
            if name.endswith(self.extrafastaExts):
                out_dict_extra_fasta[os.path.basename(name)] = self.sha256Calculator(name)
            elif name.endswith(self.genomic_ext):
                out_dict_fasta[os.path.basename(name)] = self.sha256Calculator(name)
            elif name.endswith(self.protein_ext):
                out_dict_faa[os.path.basename(name)] = self.sha256Calculator(name)
            elif name.endswith(self.allbutFasta):
                out_dict[os.path.basename(name)] = self.sha256Calculator(name)
        return (out_dict, out_dict_fasta, out_dict_faa, out_dict_extra_fasta)


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

    args = parser.parse_args()

    try:
        update_mngr = UpdateRefSeqFolder(args.new_refseq)
        update_mngr.runComparison(
            args.ftp_refseq, args.new_refseq, args.ftp_genome_dirs, args.old_genome_dirs)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise