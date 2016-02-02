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

__prog_name__ = 'verify_ani.py'
__prog_desc__ = 'Calculate intra-cluster ANI scores between representative and all genomes in cluster.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2016'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import tempfile
import shutil
from collections import defaultdict

from numpy import mean, std

from biolib.parallel import Parallel


class VerifyANI(object):
    """Calculate ANI between representative genome and all genomes in cluster."""

    def __init__(self):
        pass
    
    def _producer(self, data):
        rep_gene_file, gene_file, genome_id = data
        
        outdir = tempfile.mkdtemp()
        outfile = os.path.join(outdir, 'working_' + genome_id)
        cmd = 'ani_calculator -genome1fna %s -genome2fna %s -outfile %s -outdir %s' % (rep_gene_file, gene_file, outfile, outdir)
        os.system(cmd)
        
        with open(outfile) as f:
            f.readline()
            results = f.readline().strip().split('\t')
            gANI = 0.5*(float(results[2]) + float(results[3]))
            AF = 0.5*(float(results[4]) + float(results[5]))
        shutil.rmtree(outdir)
  
        return (genome_id, gANI, AF)

    def _consumer(self, produced_data, consumer_data):
        if consumer_data == None:
            # setup structure for consumed data
            consumer_data = []

        consumer_data.append(produced_data)

        return consumer_data

    def _progress(self, processed_items, total_items):
        return 'Processed %d of %d genomes.' % (processed_items, total_items)

    def run(self, cluster_file, genome_dir_file, output_prefix):
        output_dir = os.path.dirname(output_prefix)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # get path to all nucleotide gene file of all genomes
        gene_files = {}
        for line in open(genome_dir_file):
            line_split = line.strip().split('\t')
            
            genome_id = line_split[0]
            genome_path = line_split[1]
            genome_dir_id = os.path.basename(os.path.normpath(genome_path))
            gene_files[genome_id] = os.path.join(line_split[1], genome_dir_id + '_protein.fna')
            
        print 'Read path for %d genomes.' % len(gene_files)

        # process all clusters
        fout = open(output_prefix + '.ani.tsv', 'w')
        fout_summary = open(output_prefix + '.ani_summary.tsv', 'w')
        
        for line in open(cluster_file):
            line_split = line.strip().split('\t')
            
            rep_genome = line_split[0]
            rep_gene_file = gene_files[rep_genome]
            
            if len(line_split) == 4:
                data_items = []
                genome_ids = line_split[3].split(',')
                for genome_id in genome_ids:
                    gene_file = gene_files[genome_id]
                    data_items.append((rep_gene_file, gene_file, genome_id))
                    
                parallel = Parallel(cpus = 38)
                results = parallel.run(self._producer,
                                        self._consumer,
                                        data_items,
                                        self._progress)

                gANIs = []
                AFs = []
                for r in results:
                    genome_id, gANI, AF =  r
                    fout.write('%s\t%s\t%.3f\t%.3f\n' % (rep_genome, genome_id, gANI, AF))
                    gANIs.append(gANI)
                    AFs.append(AF)
                                            
                fout_summary.write('%s\t%.3f\t%.4f\t%.4f\t%.3f\t%.4f\t%.4f\n' % (rep_genome, 
                                                                                    mean(gANIs), 
                                                                                    std(gANIs), 
                                                                                    min(gANIs), 
                                                                                    mean(AFs), 
                                                                                    std(AFs), 
                                                                                    min(AFs)))
                fout.flush()
                fout_summary.flush()
            
        fout.close()
        fout_summary.close()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('cluster_file', help='cluster file from genometreetk aai_cluster')
    parser.add_argument('genome_dir_file', help='file indicating path to all genome directories')
    parser.add_argument('output_prefix', help='prefix for output files')
  
    args = parser.parse_args()

    try:
        p = VerifyANI()
        p.run(args.cluster_file, args.genome_dir_file, args.output_prefix)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
