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

__prog_name__ = 'validate_inter_species.py'
__prog_desc__ = 'Calculate ANI values between species assigned to the same genus.'

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
import argparse
import tempfile
import shutil
from collections import defaultdict

from numpy import mean, std

from biolib.parallel import Parallel

import dendropy


class ValidateSisterSpecies(object):
    """Calculate ANI between genomes."""

    def __init__(self):
        pass
    
    def _producer(self, data):
        genus, species_i, gids_i, species_j, gids_j = data
        
        gANIs = []
        AFs = []
        for gid_i in gids_i:
            gene_file_i = self.gene_files[gid_i]
        
            for gid_j in gids_j:
                gene_file_j = self.gene_files[gid_j]
        
                outdir = tempfile.mkdtemp()
                outfile = os.path.join(outdir, 'working_' + species_i.replace(' ', '_') + '_' + species_j.replace(' ', '_'))
                cmd = 'ani_calculator -genome1fna %s -genome2fna %s -outfile %s -outdir %s > /dev/null 2>&1' % (gene_file_i, 
                                                                                                                gene_file_j, 
                                                                                                                outfile, 
                                                                                                                outdir)
                os.system(cmd)
        
                with open(outfile) as f:
                    f.readline()
                    results = f.readline().strip().split('\t')
                    gANIs.append(0.5*(float(results[2]) + float(results[3])))
                    AFs.append(0.5*(float(results[4]) + float(results[5])))
                shutil.rmtree(outdir)
  
        return (genus, species_i, len(gids_i), species_j, len(gids_j), mean(gANIs), mean(AFs))

    def _consumer(self, produced_data, consumer_data):
        if consumer_data == None:
            consumer_data = []
            
        genus, species_i, num_species_i, species_j, num_species_j, mean_ani, mean_af = produced_data
        same_species = 'False'
        if (mean_ani >= 96.5 and mean_af >= 0.6):
            same_species = 'True'
        
        result = '%s\t%s\t%d\t%s\t%d\t%.3f\t%.3f\t%s' % (genus,
                                                            species_i,
                                                            num_species_i,
                                                            species_j,
                                                            num_species_j,
                                                            mean_ani,
                                                            mean_af,
                                                            same_species)
        consumer_data.append(result)

        return consumer_data

    def _progress(self, processed_items, total_items):
        return 'Processed %d of %d species pairs.' % (processed_items, total_items)

    def run(self, gtdb_taxonomy_file,  genome_dir_file, output_dir):
        """Calculate ANI between genomes defined as being from the same species."""
        
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # get path to all nucleotide gene files of all genomes
        self.gene_files = {}
        for line in open(genome_dir_file):
            line_split = line.strip().split('\t')
            
            gtdb_genome_id = line_split[0]
            genome_id = gtdb_genome_id.replace('GB_', '').replace('RS_', '')
            genome_path = line_split[1]
            self.gene_files[gtdb_genome_id] = os.path.join(genome_path, 'prodigal', genome_id + '_protein.fna')
            
        print('Read path for %d genomes.' % len(self.gene_files))

        # identify genomes assigned to same GTDB species
        gtdb_species = defaultdict(list)
        gtdb_genus = {}
        for line in open(gtdb_taxonomy_file):
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            taxa = [taxon.strip() for taxon in line_split[1].split(';')]
            species = gid
            genus = None
            for t in taxa:
                if t.startswith('s__') and t != 's__':
                    species = t
                elif t.startswith('g__'):
                    genus = t
                    
            if genus:
                gtdb_species[species].append(gid)
                gtdb_genus[species] = genus
                
        print('Identified %d species.' % len(gtdb_species))
 
        # process all pairs of species in the same genus
        print('Finding all pairs of species in the same genus.')
        species_list = list(gtdb_species)
        data_items = []
        for i in xrange(0, len(species_list)):
            species_i = species_list[i]
            for j in xrange(i+1, len(species_list)):
                species_j = species_list[j]

                if gtdb_genus[species_i] != gtdb_genus[species_j]:
                    continue

                data_items.append((gtdb_genus[species_i], 
                                    species_i, 
                                    gtdb_species[species_i], 
                                    species_j, 
                                    gtdb_species[species_j]))
                
        print('Identified %d pairs of species to compare.' % len(data_items))
        parallel = Parallel(cpus = 38)
        results = parallel.run(self._producer,
                                self._consumer,
                                data_items,
                                self._progress)
                                
        ani_file = os.path.join(self.output_dir, 'inter_species_ani.tsv')
        fout = open(ani_file, 'w')
        fout.write('GTDB Genus\tGTDB Species A\tNo. Genomes\tGTDB Species B\tNo. Genomes\tMean ANI\tMean AF\tSame species\n')
        for r in results:
            fout.write(r + '\n')
        fout.close()
                                        
if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_taxonomy_file', help='file specifying GTDB taxonomy of all genomes')
    parser.add_argument('genome_dir_file', help='file indicating path to all genome directories')
    parser.add_argument('output_dir', help='output directory')
  
    args = parser.parse_args()

    try:
        p = ValidateSisterSpecies()
        p.run(args.gtdb_taxonomy_file,
                args.genome_dir_file, 
                args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
