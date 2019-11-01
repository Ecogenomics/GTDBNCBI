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

__prog_name__ = 'update_genbank_from_ftp_files.py'
__prog_desc__ = ('Check which genomes are not present in Refseq but present in Genbank.' +
                 'If the genome is in Genbank, it is copied either from the ftp website for a new genome' +
                 ' or from the previous genbank folder if it\'s an existing genome')

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2016'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

import os
import sys
import shutil
import glob
import gzip
import hashlib
import argparse
import tempfile
from datetime import datetime
import multiprocessing as mp


class UpdateGenbankFolder(object):

    def __init__(self, new_genbank_folder, cpus):
        self.domains = ["archaea", "bacteria"]
        self.genome_domain_dict = {}

        self.genomic_ext = "_genomic.fna.gz"
        self.protein_ext = "_protein.faa.gz"
        self.cds_ext = "_cds_from_genomic.fna.gz"
        self.rna_ext = "_rna_from_genomic.fna.gz"

        self.fastaExts = (self.genomic_ext, self.protein_ext)
        self.extrafastaExts = (self.cds_ext, self.rna_ext)

        self.extensions = ("_feature_table.txt.gz", "_genomic.gbff.gz",
                           "_genomic.gff.gz", "_protein.gpff.gz", "_wgsmaster.gbff.gz")
        self.reports = ("_assembly_report.txt",
                        "_assembly_stats.txt", "_hashes.txt")
        self.allExts = self.fastaExts + self.extensions + self.reports
        self.allbutFasta = self.extensions + self.reports

        self.threads = cpus

        self.log = open(os.path.join(new_genbank_folder,
                                     "extra_gbk_report_gcf.log"), "w")
        self.genomes_to_review = open(os.path.join(
            new_genbank_folder, "gcaid_to_review.log"), "w", 1)
        self.select_gca = open(os.path.join(
            new_genbank_folder, "gca_selection.log"), "w")

    def runComparison(self, ftp_genbank, new_genbank, ftp_genbank_genome_dirs, old_genbank_genome_dirs, new_refseq_genome_dirs, gbk_arc_assembly, gbk_bac_assembly):
        '''
        runComparison function is walking across all directories recursively
        only folder containing latest_assembly_versions but not containing _assembly_structure
        are of interest
        '''

        # old_dict lists all records from the previous gtdb update
        with open(old_genbank_genome_dirs, 'r') as old_file:
            old_dict = {old_line.split("\t")[0]: old_line.split("\t")[1].strip()
                        for old_line in old_file}
        print("{}:old_dict loaded".format(str(datetime.now())))

        # This is a one of for the transition from the old directory structure
        # to the new one
        with open(old_genbank_genome_dirs, 'r') as old_file:
            old_dict_domain = {old_line.split("\t")[0]: old_line.split("\t")[1].strip().split("/")[8]
                               for old_line in old_file}
        print("{}:old_dict_domain loaded".format(str(datetime.now())))

        listGCA = self.parseAssemblySummary(
            gbk_arc_assembly, gbk_bac_assembly, new_refseq_genome_dirs)
        # new dict lists all records from FTP which are in new_list
        new_dict = {}
        with open(ftp_genbank_genome_dirs, 'r') as new_genome_dirs_file:
            for new_line in new_genome_dirs_file:
                new_line_split = new_line.split("\t")
                if new_line_split[0].startswith("GCA") and new_line_split[0] in listGCA:
                    new_dict[new_line_split[0]] = new_line_split[1].strip()
        print("{}:new_dict loaded".format(str(datetime.now())))

        # new genomes in FTP
        added_dict = {added_key: new_dict[added_key] for added_key in list(
            set(new_dict.keys()) - set(old_dict.keys()))}
        print "{0} genomes to add".format(len(added_dict))
        sub_added_dict = {k: added_dict[k] for k in added_dict.keys()[0:10]}
        self.addGenomes(added_dict, ftp_genbank, new_genbank)

        # delete genomes from the Database
        removed_dict = {removed_key: old_dict[removed_key] for removed_key in list(
            set(old_dict.keys()) - set(new_dict.keys()))}
        print "{0} genomes to remove".format(len(removed_dict))
        sub_removed_dict = {k: removed_dict[k]
                            for k in removed_dict.keys()[0:10]}
        self.removeGenomes(removed_dict, old_dict_domain)

        print("{}:Generating intersection list.....".format(str(datetime.now())))
        intersect_list = list(
            set(old_dict.keys()).intersection(set(new_dict.keys())))
        print("{}:Intersection list:{} genomes".format(
            str(datetime.now()), len(intersect_list)))
        self.compareGenomes(
            intersect_list, old_dict, new_dict, ftp_genbank, new_genbank, self.threads)
        self.select_gca.close()
        self.log.close()
        self.genomes_to_review.close()

    def addGenomes(self, added_dict, ftp_genbank, new_genbank):
        '''
        addGenomes function insert new genomes in the GTDB database. New genomes are present in the FTP folder
        but not in the previous version of GTDB.

        :TODO: Check if the new genome is a new version of an existing genome. in that case we overwrite the previous one
        and keep the same database id
        This will cause a conflict with the removeGenomes function.


        :param added_dict: dictionary of genomes to be added (genome_id:path to genome)
        :param ftp_genbank: base directory leading the the FTP repository for refseq
        :param new_genbank:base directory leading the new repository for refseq
        '''

        processedItems = 0
        for gcf_record in added_dict:
            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) genome to add.' % (
                processedItems, len(added_dict), float(processedItems) * 100 / len(added_dict))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            target_dir = os.path.join(
                new_genbank, added_dict[gcf_record].replace(ftp_genbank, ''))
            shutil.copytree(added_dict[
                            gcf_record], target_dir, ignore=shutil.ignore_patterns("*_assembly_structure"))
            self.log.write(
                "{0}\t{1}\tnew\n".format(self.genome_domain_dict.get(gcf_record).upper(), gcf_record))
            for compressed_file in glob.glob(target_dir + "/*.gz"):
                if os.path.isdir(compressed_file) is False:
                    inF = gzip.open(compressed_file, 'rb')
                    outF = open(
                        self.rreplace(compressed_file, ".gz", "", 1), 'wb')
                    outF.write(inF.read())
                    inF.close()
                    outF.close()
                    os.remove(compressed_file)
        sys.stdout.write('\n')

    def removeGenomes(self, removed_dict, old_dict_domain):
        '''
        removeGenomes function removes all outdated genomes from the gtdb database
        In addition it tracks the lists(name and owner) that have been modified while deleting those genomes

        :param removed_dict: dictionary of genomes to delete
        '''

        for gca_record in removed_dict:
            self.log.write(
                "{0}\t{1}\tremoved\n".format(old_dict_domain.get(gca_record).upper(), gca_record))

    def compareGenomes(self, intersect_list, old_dict, new_dict, ftp_genbank, new_genbank, threads):
        '''
        compare the genomes existing in both folders ( FTP folder and previous gtdb update).

        :param intersect_list:
        :param old_dict:
        :param new_dict:
        '''

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for gca_record in intersect_list:
            gtdb_dir = old_dict.get(gca_record)
            ftp_dir = new_dict.get(gca_record)
            target_dir = os.path.join(
                new_genbank, ftp_dir.replace(ftp_genbank, ''))
            workerQueue.put((gtdb_dir, ftp_dir, target_dir, gca_record))

        for _ in range(threads):
            workerQueue.put((None, None, None, None))

        try:
            workerProc = [mp.Process(target=self.__workerThread, args=(
                workerQueue, writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target=self.__writerThread, args=(
                len(intersect_list), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()

            writeProc.terminate()

    def __workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            gtdb_dir, ftp_dir, target_dir, gca_record = queueIn.get(
                block=True, timeout=None)
            if gca_record is None:
                break

            status_gca = self.readmd5Checksum(
                gtdb_dir, ftp_dir, target_dir, gca_record)
            queueOut.put(status_gca)

    def __writerThread(self, numgenometoprocess, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a is None:
                break
            self.log.write(a)

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) genome pairs.' % (
                processedItems, numgenometoprocess, float(processedItems) * 100 / numgenometoprocess)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def readmd5Checksum(self, gtdb_dir, ftp_dir, target_dir, gcf_record):
        '''
        Compare the checksum of the file listed in the checksums.txt
        '''
        pathftpmd5 = os.path.join(ftp_dir, "md5checksums.txt")
        pathgtdbmd5 = os.path.join(gtdb_dir, "md5checksums.txt")
        target_pathnewmd5 = os.path.join(target_dir, "md5checksums.txt")
        status = []

        tmp_ftp_dir = tempfile.mkdtemp()
        tmp_target = os.path.join(tmp_ftp_dir, os.path.basename(target_dir))
        shutil.copytree(ftp_dir, tmp_target, symlinks=True,
                        ignore=shutil.ignore_patterns("*_assembly_structure"))
        for compressed_file in glob.glob(tmp_target + "/*.gz"):
            if os.path.isdir(compressed_file) is False:
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

        ftpdict, ftpdict_fasta, ftpdict_faa, ftpdict_extra_fasta = self.parse_checksum(
            tmp_target)
        gtdbdict, gtdbdict_fasta, gtdbdict_faa, gtdbdict_extra_fasta = self.parse_checksum(
            gtdb_dir)

        # if the genomic.fna.gz or the protein.faa.gz are missing, we set this
        # record as incomplete
        if len(list(set(ftpdict_fasta.keys()).symmetric_difference(set(gtdbdict_fasta.keys())))) > 0:
            self.genomes_to_review.write("tmp_target:{}\nftpdict_fasta.keys():{}\ngtdb_dir:{}\ngtdbdict_fasta.keys():{}\n\n".format(tmp_target,
                                                                                                                                    ftpdict_fasta.keys(),
                                                                                                                                    gtdb_dir,
                                                                                                                                    gtdbdict_fasta.keys()))
            status.append("incomplete")
            shutil.copytree(tmp_target, target_dir, symlinks=True,
                            ignore=shutil.ignore_patterns("*_assembly_structure"))
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
                    elif len(target_files) == 0 and len(ftp_files) == 0 and report == '_hashes.txt':
                        status.append("old_folder_dir")
                    elif len(target_files) == 0 and len(ftp_files) == 1 and report == '_hashes.txt':
                        shutil.copy2(ftp_dir[0], ftp_dir[0].replace(
                            ftp_dir.target_dir))
                        status.append("new_hashes")
                    else:
                        print "########"
                        print target_dir
                        print target_files
                        print ftp_dir
                        print ftp_files
                        print "IT SHOULDN'T HAPPEN"
                        print "########"
                        sys.exit()
        status_gca = "{0}\t{1}\t{2}\n".format(self.genome_domain_dict.get(
            gcf_record).upper(), gcf_record, ';'.join([x for x in set(status)]))
        shutil.rmtree(tmp_ftp_dir)
        return status_gca

# Tools

    def parseAssemblySummary(self, gbk_arc_assembly, gbk_bac_assembly, new_refseq_genome_dirs):
        listGCA = []
        dictGCF = self._populateGenomesDict(new_refseq_genome_dirs)
        print "parsing of dictionary is done"

        for domain, assemblyfile in [('archaea', gbk_arc_assembly), ('bacteria', gbk_bac_assembly)]:
            num_lines = sum(1 for line in open(assemblyfile))
            processedItems = 0
            with open(assemblyfile, "r") as sumf:

                # we discard the first line
                sumf.readline()
                for line in sumf:
                    processedItems += 1
                    statusStr = 'Finished processing %d of %d (%.2f%%) gca records.' % (
                        processedItems, num_lines, float(processedItems) * 100 / num_lines)
                    sys.stdout.write('%s\r' % statusStr)
                    sys.stdout.flush()
                    split_line = line.split("\t")
                    gcf_access = 'G' + split_line[17][4:13]
                    full_gca_access = split_line[0]
                    latest = split_line[10]
                    surveillance_info = split_line[20]
                    if 'surveillance' in surveillance_info:
                        continue
                    elif latest == "latest":
                        if not gcf_access.startswith("GCF"):
                            formatted_gcaid = 'G' + full_gca_access[4:13]
                            if formatted_gcaid in dictGCF:
                                self.select_gca.write("{0} skipped because {1} in Refseq (although {0} has no GCF)\n".format(
                                    full_gca_access, gcf_access))
                                continue
                            else:
                                listGCA.append(full_gca_access)
                        else:
                            # if the Refseq folder is empty, we copy the genbank
                            # folder
                            if gcf_access in dictGCF:
                                protein_files = glob.glob(
                                    os.path.join(dictGCF.get(gcf_access), "*_protein.faa"))
                                if len(protein_files) == 0:
                                    self.select_gca.write(
                                        "{0} associated with {1} : {1} missed files in FTP folder\n".format(full_gca_access, gcf_access))
                                    listGCA.append(full_gca_access)
                            else:
                                self.select_gca.write(
                                    "{0} associated with {1} : {1} not present in FTP folder\n".format(full_gca_access, gcf_access))
                                listGCA.append(full_gca_access)
                    self.genome_domain_dict[full_gca_access] = domain
        sys.stdout.write('\n')
        return listGCA

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

    def _populateGenomesDict(self, genome_dirs_file):
        temp_dict = {}
        with open(genome_dirs_file, "r") as list_dirs:
            for line in list_dirs:
                temp_dict['G' + line.split("\t")[0][4:13]] = line.split(
                    "\t")[1].rstrip()
        return temp_dict

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

    def parse_checksum(self, pathtodir):
        '''
        parse_checksum function parses the md5 checksum file.
        It returns 2 dictionaries {file:size} : one for the fna and faa files, one for the genbank files
        :param md5File:
        '''

        out_dict, out_dict_fasta, out_dict_faa, out_dict_extra_fasta = {}, {}, {}, {}

        for name in glob.glob(os.path.join(pathtodir, '*')):
            if name.endswith(self.extrafastaExts):
                out_dict_extra_fasta[os.path.basename(
                    name)] = self.sha256Calculator(name)
                os.chmod(name, 0o664)
            elif name.endswith(self.genomic_ext):
                out_dict_fasta[os.path.basename(
                    name)] = self.sha256Calculator(name)
                os.chmod(name, 0o664)
            elif name.endswith(self.protein_ext):
                out_dict_faa[os.path.basename(
                    name)] = self.sha256Calculator(name)
                os.chmod(name, 0o664)
            elif name.endswith(self.allbutFasta):
                out_dict[os.path.basename(name)] = self.sha256Calculator(name)
                os.chmod(name, 0o664)
        return (out_dict, out_dict_fasta, out_dict_faa, out_dict_extra_fasta)


if __name__ == "__main__":
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ftp_genbank_directory', dest="ftp_genbank", required=True,
                        help='base directory leading the the FTP repository for genbank')
    parser.add_argument('--new_genbank_directory', dest="new_genbank",
                        required=True, help='base directory leading the new repository for genbank')
    parser.add_argument('--ftp_genbank_genome_dirs_file', dest="ftp_genbank_genome_dirs", required=True,
                        help='metadata file listing all directories for the FTP folder (generated by ncbi_genome_dirs.py)')
    parser.add_argument('--old_genbank_genome_dirs_file', dest="old_genbank_genome_dirs", required=True,
                        help='metadata file listing all directories from the previous NCBI update date  (generated by genome_dirs.py)')
    parser.add_argument('--new_refseq_genome_dirs_file', dest="new_refseq_genome_dirs", required=True,
                        help='metadata file listing all directories from the previous NCBI update date  (generated by genome_dirs.py)')
    parser.add_argument('--arc_assembly_summary', required=True,
                        help='Genbank metadata file downloaded from ncbi.')
    parser.add_argument('--bac_assembly_summary', required=True,
                        help='Genbank metadata file downloaded from ncbi.')
    parser.add_argument('--cpus', type=int, default=1,
                        help='Number of cpus')
    args = parser.parse_args()

    try:
        update_manager = UpdateGenbankFolder(args.new_genbank, args.cpus)
        update_manager.runComparison(
            args.ftp_genbank, args.new_genbank, args.ftp_genbank_genome_dirs,
            args.old_genbank_genome_dirs, args.new_refseq_genome_dirs,
            args.arc_assembly_summary, args.bac_assembly_summary)

    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
