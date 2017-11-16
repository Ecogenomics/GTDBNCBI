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

__prog_name__ = 'validate_species.py'
__prog_desc__ = 'Calculate ANI values between genomes classified as being from the same species.'

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


class ValidateANI(object):
    """Calculate ANI between genomes defined as being from the same species."""

    def __init__(self):
        pass
    
    def _producer(self, data):
        genome_id_a, genome_id_b, gene_file_a, gene_file_b = data
        
        outdir = tempfile.mkdtemp()
        outfile = os.path.join(outdir, 'working_' + genome_id_a + '_' + genome_id_b)
        cmd = 'ani_calculator -genome1fna %s -genome2fna %s -outfile %s -outdir %s' % (gene_file_a, gene_file_b, outfile, outdir)
        os.system(cmd)
        
        with open(outfile) as f:
            f.readline()
            results = f.readline().strip().split('\t')
            gANI = 0.5*(float(results[2]) + float(results[3]))
            AF = 0.5*(float(results[4]) + float(results[5]))
        shutil.rmtree(outdir)
  
        return (genome_id_a, genome_id_b, gANI, AF)

    def _consumer(self, produced_data, consumer_data):
        if consumer_data == None:
            # setup structure for consumed data
            consumer_data = []

        consumer_data.append(produced_data)

        return consumer_data

    def _progress(self, processed_items, total_items):
        return 'Processed %d of %d genome pairs.' % (processed_items, total_items)

    def run(self, gtdb_taxonomy_file, gtdb_metadata_file, genome_dir_file, output_dir):
        """Calculate ANI between genomes defined as being from the same species."""
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # get path to all nucleotide gene files of all genomes
        gene_files = {}
        for line in open(genome_dir_file):
            line_split = line.strip().split('\t')
            
            gtdb_genome_id = line_split[0]
            genome_id = gtdb_genome_id.replace('GB_', '').replace('RS_', '')
            genome_path = line_split[1]
            gene_files[gtdb_genome_id] = os.path.join(genome_path, 'prodigal', genome_id + '_protein.fna')
            
        print('Read path for %d genomes.' % len(gene_files))
        
        # get NCBI species designations
        ncbi_species = {}
        ncbi_type_strain = set()
        lpsn_type_species = set()
        with open(gtdb_metadata_file) as f:
            header = f.readline().strip().split('\t')
            
            ncbi_taxonomy_index = header.index('ncbi_taxonomy')
            ncbi_type_strain_index = header.index('ncbi_type_strain')
            lpsn_type_species_index = header.index('lpsn_species')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[0]
                taxonomy = [taxon.strip() for taxon in line_split[ncbi_taxonomy_index].split(';')]
                if len(taxonomy) == 7:
                    sp = taxonomy[6]
                    if sp != 's__':
                        ncbi_species[gid] = sp
                        
                        if line_split[ncbi_type_strain_index] == 'yes':
                            ncbi_type_strain.add(sp)
                            
                        if line_split[lpsn_type_species_index] == sp:
                            lpsn_type_species.add(sp)

        # identify genomes assigned to same GTDB species
        gtdb_species = defaultdict(list)
        for line in open(gtdb_taxonomy_file):
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            taxonomy = [taxon.strip() for taxon in line_split[1].split(';')]
            sp = taxonomy[6]
            if sp != 's__':
                gtdb_species[sp].append(gid)
                
        print('Identified %d species.' % len(gtdb_species))
        
        # count species with multiple genomes
        multi_species_count = 0
        for sp, genome_ids in gtdb_species.items():
            if len(genome_ids) >= 2:
                multi_species_count += 1

        print('Identified %d GTDB species with multiple genomes.' % multi_species_count)

        # process all clusters
        ani_file = os.path.join(output_dir, 'ani.tsv')
        fout = open(ani_file, 'w')
        fout.write('GTDB Species\tNo. Genomes\tGenome ID A\tGenome ID B\tgANI\tAF\tNCBI Species A\tNCBI Species B\tAccepted Species\n')
        
        summary_file = os.path.join(output_dir, 'species_ani.tsv')
        fout_summary = open(summary_file, 'w')
        fout_summary.write('GTDB Species\tNo. Genomes\tMean gANI\tStd gANI\tMin gANI\tMean AF\tStd AF\tMin AF')
        fout_summary.write('\tPass ANI Species Test\tNCBI Species\tNo. NCBI Species\tNCBI Type Strains')
        fout_summary.write('\tLPSN Type Species at NCBI\tNo. LPSN Type Species at NCBI\n')
        
        for sp, genome_ids in gtdb_species.items():
            if len(genome_ids) <= 1:
                continue
                
            data_items = []
            for a in xrange(0, len(genome_ids)):
                for b in xrange(a+1, len(genome_ids)):
                    gene_file_a = gene_files[genome_ids[a]]
                    gene_file_b = gene_files[genome_ids[b]]
                    data_items.append((genome_ids[a], genome_ids[b], gene_file_a, gene_file_b))
                
            parallel = Parallel(cpus = 38)
            results = parallel.run(self._producer,
                                    self._consumer,
                                    data_items,
                                    self._progress)

            gANIs = []
            AFs = []
            all_ncbi_sp = set()
            all_ncbi_type_strains = set()
            all_lpsn_type_species = set()
            for r in results:
                genome_id_a, genome_id_b, gANI, AF =  r
                ncbi_sp_a = ncbi_species.get(genome_id_a, 'unassigned')
                ncbi_sp_b = ncbi_species.get(genome_id_b, 'unassigned')
                all_ncbi_sp.add(ncbi_sp_a)
                all_ncbi_sp.add(ncbi_sp_b)
                
                if ncbi_sp_a in ncbi_type_strain:
                    all_ncbi_type_strains.add(ncbi_sp_a)
                if ncbi_sp_b in ncbi_type_strain:
                    all_ncbi_type_strains.add(ncbi_sp_b)
                
                if ncbi_sp_a in lpsn_type_species:
                    all_lpsn_type_species.add(ncbi_sp_a)
                if ncbi_sp_b in lpsn_type_species:
                    all_lpsn_type_species.add(ncbi_sp_b)
                
                
                accepted_species = 'invalid'
                if (gANI >= 96.5 and AF >= 0.6) or (ncbi_sp_a != 'unassigned' and (ncbi_sp_a == ncbi_sp_b)):
                    accepted_species = 'valid'
                fout.write('%s\t%d\t%s\t%s\t%.3f\t%.3f\t%s\t%s\t%s\n' % (sp, 
                                                                            len(genome_ids), 
                                                                            genome_id_a, 
                                                                            genome_id_b, 
                                                                            gANI, 
                                                                            AF,
                                                                            ncbi_sp_a,
                                                                            ncbi_sp_b,
                                                                            accepted_species))
                gANIs.append(gANI)
                AFs.append(AF)
                
            print('%s\t%d\t%.2f\t%.2f' % (sp, len(genome_ids), mean(gANIs), min(gANIs)))
            ani_species = 'False'
            if (min(gANIs) >= 96.5 and min(AFs) >= 0.6):
                ani_species = 'True'
            fout_summary.write('%s\t%d\t%.3f\t%.4f\t%.4f\t%.3f\t%.4f\t%.4f\t%s\t%s\t%d\t%s\t%s\t%d\n' % (
                                    sp,
                                    len(genome_ids),
                                    mean(gANIs), 
                                    std(gANIs), 
                                    min(gANIs), 
                                    mean(AFs), 
                                    std(AFs), 
                                    min(AFs),
                                    ani_species,
                                    ','.join(all_ncbi_sp),
                                    len(all_ncbi_sp) - int('unassigned' in all_ncbi_sp),
                                    ','.join(all_ncbi_type_strains),
                                    ','.join(all_lpsn_type_species),
                                    len(all_lpsn_type_species)
                                    ))
            fout.flush()
            fout_summary.flush()
            
        fout.close()
        fout_summary.close()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtdb_taxonomy_file', help='file specifying GTDB taxonomy of all genomes')
    parser.add_argument('gtdb_metadata_file', help='file specifying GTDB metadata (TSV)')
    parser.add_argument('genome_dir_file', help='file indicating path to all genome directories')
    parser.add_argument('output_dir', help='output directory')
  
    args = parser.parse_args()

    try:
        p = ValidateANI()
        p.run(args.gtdb_taxonomy_file, 
                args.gtdb_metadata_file, 
                args.genome_dir_file, 
                args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
