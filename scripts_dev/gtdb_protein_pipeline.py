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

__prog_name__ = 'gtdb_protein_pipeline.py'
__prog_desc__ = 'Run genomes through Prodigal, Pfam and TIGRfam annotation.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2018'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import shutil
import ntpath
import argparse
import logging
from collections import defaultdict
import multiprocessing as mp

from biolib.logger import logger_setup
from biolib.common import remove_extension
from biolib.external.execute import check_dependencies
from biolib.external.prodigal import Prodigal
from biolib.checksum import sha256


class Pipeline(object):
    """TBD"""

    def __init__(self, tmp_dir, output_dir, cpus):
        """Initialization."""
        
        self.tmp_dir = tmp_dir
        self.output_dir = output_dir
        self.cpus = cpus
        
        check_dependencies(['prodigal', 'hmmsearch', 'pfam_search.pl', 'genometk'])

        self.tigrfam_hmms = '/srv/whitlam/bio/db/tigrfam/15.0/TIGRFAMs_15.0_HMM/tigrfam.hmm'
        self.tigrfam_ext = '_tigrfam.tsv'
        
        self.pfam_hmm_dir = '/srv/db/pfam/27/'
        self.pfam_ext = '_pfam.tsv'
        
        self.protein_file_ext = '_protein.faa'
        
        logger_setup(output_dir, "gtdb_protein_pipeline.log", "gtdb_protein_pipeline", __version__, False)
        self.logger = logging.getLogger('timestamp')
    
    def _parse_translation_table(self, gene_feature_file):
        """Parse translation table from NCBI GFF file."""

        ncbi_transl_table = None
        if os.path.exists(gene_feature_file):
            for line in open(gene_feature_file):
                if line[0] == '#':
                    continue
                    
                if 'transl_table=' in line:
                    trans_num = line[line.rfind('=')+1:].strip()
                    gene_feature_file = int(trans_num)
                    break
                    
        return gene_feature_file
        
    def _run_prodigal(self, genome_paths):
        """Run Prodigal on genomes."""

        # get genome path and translation table for each file
        self.logger.info('Determining genomic file and translation table for each of the %d genomes.' % len(genome_paths))
        genome_files = []
        translation_table = {}
        for gid, gpath in genome_paths.items():
            assembly_id = os.path.basename(os.path.normpath(gpath))
            canonical_gid = assembly_id[0:assembly_id.find('_', 4)]
            
            genome_file = os.path.join(gpath, assembly_id + '_genomic.fna')
            if os.path.exists(genome_file):
                if os.stat(genome_file).st_size == 0:
                    self.logger.warning('Genomic file appears to be empty: %s' % genome_file)
                    continue
                
                genome_files.append(genome_file)
            else:
                self.logger.warning('Genomic file appears to be missing: %s' % genome_file)
                    
            gff_file = os.path.join(gpath, assembly_id + '_genomic.gff')
            if os.path.exists(gff_file):
                if os.stat(gff_file).st_size == 0:
                    self.logger.warning('GFF appears to be empty: %s' % gff_file)
                    continue

                tt = self._parse_translation_table(gff_file)
                if tt:
                    translation_table[canonical_gid] = tt
                else:
                    translation_table[canonical_gid] = None
                    self.logger.warning('Unable to determine translation table for: %s' % gff_file)
                    sys.exit(-1)
            else:
                self.logger.warning('GFF appears to be missing: %s' % gff_file)
                sys.exit(-1)
        
        # run Prodigal on each genome
        self.logger.info('Running Prodigal on %d genomes.' % len(genome_paths))
        prodigal = Prodigal(cpus=self.cpus)
        summary_stats = prodigal.run(genome_files, 
                                    translation_table=translation_table, 
                                    output_dir=self.tmp_dir)

        # move results into individual genome directories
        self.logger.info('Moving files and calculating checksums.')
        for genome_file in genome_files:
            genome_path, genome_id = ntpath.split(genome_file)
            genome_id = remove_extension(genome_id)
            canonical_gid = genome_id[0:genome_id.find('_', 4)]
            
            aa_gene_file = os.path.join(self.tmp_dir, genome_id + '_genes.faa')
            nt_gene_file = os.path.join(self.tmp_dir, genome_id + '_genes.fna')
            gff_file = os.path.join(self.tmp_dir, genome_id + '.gff')

            genome_root = genome_id[0:genome_id.find('_', 4)]
            prodigal_path = os.path.join(genome_path, 'prodigal')
            if not os.path.exists(prodigal_path):
                os.makedirs(prodigal_path)
            new_aa_gene_file = os.path.join(prodigal_path, genome_root + '_protein.faa')
            new_nt_gene_file = os.path.join(prodigal_path, genome_root + '_protein.fna')
            new_gff_file = os.path.join(prodigal_path, genome_root + '_protein.gff')

            os.system('mv %s %s' % (aa_gene_file, new_aa_gene_file))
            os.system('mv %s %s' % (nt_gene_file, new_nt_gene_file))
            os.system('mv %s %s' % (gff_file, new_gff_file))

            # save translation table information
            translation_table_file = os.path.join(prodigal_path, 'prodigal_translation_table.tsv')
            fout = open(translation_table_file, 'w')
            if translation_table[canonical_gid]:
                fout.write('%s\t%d\t%s\n' % ('best_translation_table', 
                                                summary_stats[genome_id].best_translation_table,
                                                'used table specified by NCBI'))
            else:
                fout.write('%s\t%d\n' % ('best_translation_table', summary_stats[genome_id].best_translation_table))
                fout.write('%s\t%.2f\n' % ('coding_density_4', summary_stats[genome_id].coding_density_4 * 100))
                fout.write('%s\t%.2f\n' % ('coding_density_11', summary_stats[genome_id].coding_density_11 * 100))
            fout.close()

            checksum = sha256(new_aa_gene_file)
            fout = open(new_aa_gene_file + '.sha256', 'w')
            fout.write(checksum)
            fout.close()
            
    def _tigr_top_hit(self, tigrfam_file, tigrfam_tophit_file):
        """Identify top TIGRfam hits."""

        tophits = {}
        for line in open(tigrfam_file):
            if line[0] == '#' or line[0] == '[':
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[3]
            evalue = float(line_split[4])
            bitscore = float(line_split[5])
            if gene_id in tophits:
                if bitscore > tophits[gene_id][2]:
                    tophits[gene_id] = (hmm_id, evalue, bitscore)
            else:
                tophits[gene_id] = (hmm_id, evalue, bitscore)

        fout = open(tigrfam_tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, stats in tophits.iteritems():
            hit_str = ','.join(map(str, stats))
            fout.write('%s\t%s\n' % (gene_id, hit_str))
        fout.close()

        # calculate checksum
        checksum = sha256(tigrfam_tophit_file)
        fout = open(tigrfam_tophit_file + '.sha256', 'w')
        fout.write(checksum)
        fout.close()
            
    def __tigr_worker(self, queue_in, queue_out):
        """Process each data item in parallel."""
        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file == None:
                break

            assembly_dir, filename = os.path.split(gene_file)

            output_hit_file = os.path.join(assembly_dir, filename.replace(self.protein_file_ext, '_tigrfam.tsv'))
            hmmsearch_out = os.path.join(assembly_dir, filename.replace(self.protein_file_ext, '_tigrfam.out'))
            cmd = 'hmmsearch -o %s --tblout %s --noali --notextw --cut_nc --cpu 1 %s %s' % (hmmsearch_out, output_hit_file, self.tigrfam_hmms, gene_file)
            os.system(cmd)

            # calculate checksum
            checksum = sha256(output_hit_file)
            fout = open(output_hit_file + '.sha256', 'w')
            fout.write(checksum)
            fout.close()
            
            # determine top hits
            tigrfam_tophit_file = os.path.join(assembly_dir, filename.replace(self.protein_file_ext, '_tigrfam_tophit.tsv'))
            self._tigr_top_hit(output_hit_file, tigrfam_tophit_file)

            # allow results to be processed or written to file
            queue_out.put(gene_file)
            
    def _pfam_top_hit(self, pfam_file, pfam_tophit_file):
        """Identify top Pfam hits."""
        
        tophits = defaultdict(dict)
        for line in open(pfam_file):
            if line[0] == '#' or not line.strip():
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[5]
            evalue = float(line_split[12])
            bitscore = float(line_split[11])
            if gene_id in tophits:
                if hmm_id in tophits[gene_id]:
                    if bitscore > tophits[gene_id][hmm_id][1]:
                        tophits[gene_id][hmm_id] = (evalue, bitscore)
                else:
                    tophits[gene_id][hmm_id] = (evalue, bitscore)
            else:
                tophits[gene_id][hmm_id] = (evalue, bitscore)

        fout = open(pfam_tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, hits in tophits.iteritems():
            hit_str = []
            for hmm_id, stats in hits.iteritems():
                hit_str.append(hmm_id + ',' + ','.join(map(str, stats)))
            fout.write('%s\t%s\n' % (gene_id, ';'.join(hit_str)))
        fout.close()

        # calculate checksum
        checksum = sha256(pfam_tophit_file)
        fout = open(pfam_tophit_file + '.sha256', 'w')
        fout.write(checksum)
        fout.close()
            
    def __pfam_worker(self, queue_in, queue_out):
        """Process each data item in parallel."""
        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file == None:
                break

            assembly_dir, filename = os.path.split(gene_file)

            output_hit_file = os.path.join(assembly_dir, filename.replace(self.protein_file_ext, '_pfam.tsv'))
            cmd = 'pfam_search.pl -outfile %s -cpu 1 -fasta %s -dir %s' % (output_hit_file, gene_file, self.pfam_hmm_dir)
            os.system(cmd)

            # calculate checksum
            checksum = sha256(output_hit_file)
            fout = open(output_hit_file + '.sha256', 'w')
            fout.write(checksum)
            fout.close()
            
            # determine top hits
            pfam_tophit_file = os.path.join(assembly_dir, filename.replace(self.protein_file_ext, '_pfam_tophit.tsv'))
            self._pfam_top_hit(output_hit_file, pfam_tophit_file)

            # allow results to be processed or written to file
            queue_out.put(gene_file)

    def __progress(self, num_items, queue_out):
        """Store or write results of worker threads in a single thread."""
        processed_items = 0
        while True:
            a = queue_out.get(block=True, timeout=None)
            if a == None:
                break

            processed_items += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processed_items, num_items, float(processed_items)*100/num_items)
            sys.stdout.write('%s\r' % statusStr)
        sys.stdout.flush()

        sys.stdout.write('\n')
        
    def _run_annotation(self, genome_paths, worker_func, db_str):
        """Determine TIGRfam or Pfam annotations."""

        gene_files = []
        for gid, gpath in genome_paths.items():
            prodigal_dir = os.path.join(gpath, 'prodigal')
            assembly_id = os.path.basename(os.path.normpath(gpath))
            canonical_gid = assembly_id[0:assembly_id.find('_', 4)]
            gene_file = os.path.join(prodigal_dir, canonical_gid + self.protein_file_ext)
            if os.path.exists(gene_file):
                if os.stat(gene_file).st_size == 0:
                    self.logger.warning('Protein file appears to be empty: %s' % gene_file)
                else:
                    gene_files.append(gene_file)
            else:
                self.logger.warning('Protein file appears to be missing: %s' % gene_file)
        
        # populate worker queue with data to process
        self.logger.info('Determining %s annotations for %d genomes using %d CPUs.' % (db_str, len(gene_files), self.cpus))
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in gene_files:
            workerQueue.put(f)

        for _ in range(self.cpus):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target = worker_func, args = (workerQueue, writerQueue)) for _ in range(self.cpus)]
            writeProc = mp.Process(target = self.__progress, args = (len(gene_files), writerQueue))

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
            
    def _protein_metadata(self, genome_paths):
        """Calculate metadata for proteins."""
        
        self.logger.info('Calculating protein metadata for %d genomes.' % len(genome_paths))
        for gid, gpath in genome_paths.items():
            assembly_id = os.path.basename(os.path.normpath(gpath))
            canonical_gid = assembly_id[0:assembly_id.find('_', 4)]
            genome_file = os.path.join(gpath, assembly_id + '_genomic.fna')
            gff_file = os.path.join(gpath, 'prodigal', canonical_gid + '_protein.gff')

            # clean up old log files
            log_file = os.path.join(gpath, 'genometk.log')
            if os.path.exists(log_file):
                os.remove(log_file)

            # calculate metadata
            os.system('genometk gene --silent %s %s %s' % (genome_file, gff_file, gpath))

    def run(self, genome_list_file, gtdb_genome_path_file):
        """Process genomes."""
        
        # get genomes to process
        genomes_to_process = set()
        for line in open(genome_list_file):
            if line[0] != '#':
                genomes_to_process.add(line.strip().split('\t')[0])
        self.logger.info('Identified %d genomes to process.' % len(genomes_to_process))
            
        # get path of genomes to process
        genome_paths = {}
        for line in open(gtdb_genome_path_file):
            line_split = line.strip().split('\t')
            gid = line_split[0]
            if gid in genomes_to_process:
                gpath = line_split[1]
                genome_paths[gid] = gpath
                
                #***prodigal_dir = os.path.join(gpath, 'prodigal')
                #***if os.path.exists(prodigal_dir):
                #***    shutil.rmtree(prodigal_dir)

        # run Prodigal
        #***self._run_prodigal(genome_paths)
        
        # run TIGRfam and Pfam
        #***self._run_annotation(genome_paths, self.__tigr_worker, 'TIGRfam')
        #***self._run_annotation(genome_paths, self.__pfam_worker, 'Pfam')
        
        # calculate protein metadata
        self._protein_metadata(genome_paths)
        
        self.logger.info('Done.')

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_list_file', help='file listing genomes to process')
    parser.add_argument('gtdb_genome_path_file', help='file indicating path to GTDB genomes')
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('--tmp_dir', help='temporary directory for storing intermediate results', default='/tmp')
    parser.add_argument('-c', '--cpus', type=int, help='number of CPUs', default=1)

    args = parser.parse_args()

    try:
        p = Pipeline(args.tmp_dir, args.output_dir, args.cpus)
        p.run(args.genome_list_file, args.gtdb_genome_path_file)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
