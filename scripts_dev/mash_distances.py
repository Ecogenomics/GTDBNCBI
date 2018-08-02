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

__prog_name__ = 'mash_distances.py'
__prog_desc__ = 'Calculate pairwise Mash distances between GTDB genomes.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2017'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import shutil
import tempfile
import argparse

class Mash(object):
    """Calculate pairwise Mash distances between GTDB genomes."""

    def __init__(self):
        """Initialization."""        
        pass

    def run(self, genome_path_file, max_dist, max_pvalue, output_dir, threads):
        """Calculate pairwise Mash distances between GTDB genomes."""

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

        # get path to genome files
        print 'Getting path to genome files.'
        genome_files = {}
        for line in open(genome_path_file):
            line_split = line.strip().split('\t')

            gtdb_id = line_split[0]
            genome_dir = line_split[1]
            genome_acc = os.path.basename(os.path.normpath(genome_dir))

            genome_file = os.path.join(genome_dir, genome_acc + '_genomic.fna')
            genome_files[gtdb_id] = genome_file


        # create file specifying all genome files
        print 'Creating file pointing to each genome file.'
        genome_list_file = os.path.join(output_dir, 'genome_files.lst')
        fout = open(genome_list_file, 'w')
        for i, (genome_id, genome_file) in enumerate(genome_files.iteritems()):
            fout.write(genome_file + '\n')
        fout.close()
            
        # calculate pairwise Mash distance between genomes
        print 'Creating Mash sketch.'
        sketch = os.path.join(output_dir, 'gtdb.msh')
        cmd = 'mash sketch -l -p %d -k 16 -s 5000 -o %s %s' % (threads, sketch, genome_list_file)
        os.system(cmd)
        
        print 'Calculating pairwise Mash distances.'
        mash_output = os.path.join(output_dir, 'gtdb_mash_dist.tmp')
        cmd = 'mash dist -p %d -d %f -v %f %s %s > %s' % (threads,
                                                            max_dist,
                                                            max_pvalue,
                                                            sketch, 
                                                            sketch, 
                                                            mash_output)
        os.system(cmd)
        
        print 'Fixing genome names in output file.'
        fout = open(os.path.join(output_dir, 'gtdb_mash_dist.tsv'), 'w')
        for line in open(mash_output):
            line_split = line.strip().split('\t')
            
            genome_a = line_split[0]
            genome_a = genome_a[genome_a.rfind('/')+1:]
            underscore = genome_a.find('_', 4)
            if underscore != -1:
                genome_a = genome_a[0:underscore]
            if genome_a.startswith('GCF_'):
                genome_a = 'RS_' + genome_a
            elif genome_a.startswith('GCA_'):
                genome_a = 'GB_' + genome_a
            
            genome_b = line_split[1]
            genome_b = genome_b[genome_b.rfind('/')+1:]
            underscore = genome_b.find('_', 4)
            if underscore != -1:
                genome_b = genome_b[0:underscore]
            if genome_b.startswith('GCF_'):
                genome_b = 'RS_' + genome_b
            elif genome_b.startswith('GCA_'):
                genome_b = 'GB_' + genome_b
            
            fout.write('%s\t%s\t%s\n' % (genome_a, genome_b, '\t'.join(line_split[2:])))
            
        fout.close()
        
        os.remove(mash_output)

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_path_file', help='file indicating directories of GTDB genomes')
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('-d', '--max_dist', type=float, help='maximum distance to report', default=0.1)
    parser.add_argument('-p', '--max_pvalue', type=float, help='maximum p-value to report', default=1e-5)
    parser.add_argument('-t', '--threads', type=int, help='number of threads', default=1)

    args = parser.parse_args()

    try:
        p = Mash()
        p.run(args.genome_path_file, 
                args.max_dist, 
                args.max_pvalue, 
                args.output_dir, 
                args.threads)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
