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

__prog_name__ = 'dereplicate_taxa.py'
__prog_desc__ = 'Sample taxa based on GTDB taxonomy, presence of key marker genes, and NCBI metadata.'

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
import ntpath
import csv
import random
import operator
from collections import defaultdict, namedtuple

from biolib.taxonomy import Taxonomy
from biolib.seq_io import read_seq

csv.field_size_limit(sys.maxsize)

   
class DereplicateTaxa(object):
    """Sample taxa."""

    def __init__(self, min_n50, max_contig_count, max_scaffold_count):
        
        self.max_scaffolds = max_scaffold_count
        self.max_contigs = max_contig_count
        self.min_n50_scaffolds = min_n50
        self.max_ambiguous = 100000
        self.max_total_gap_length = 1e6
        
    def _read_marker_count_file(self, count_file):
        """Read marker count for genomes."""
        
        count = {}
        with open(count_file) as f:
            header = f.readline().split('\t')
            comp_index = header.index('Completeness')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                genome_id = line_split[0]
                single_copy_count = sum([1 for c in line_split[2:comp_index] if int(c) == 1])
                count[genome_id] = single_copy_count
                
        return count
        
    def _read_ssu_file(self, ssu_fasta_file):
        """Read length of SSU sequences for genomes."""
        
        ssu_length = {}
        for seq_id, seq in read_seq(ssu_fasta_file):
            genome_id = seq_id.split('~')[0]
            
            ssu_length[genome_id] = len(seq) - seq.upper().count('N')
            
        return ssu_length
        
    def _read_metadata(self, metadata_file, taxonomy_file):
        """Read metadata for genomes."""
        
        gtdb_taxonomy_from_file = {}
        if taxonomy_file:
            print('Reading taxonomy from: %s' % taxonomy_file)
            for line in open(taxonomy_file):
                line_split = line.strip().split('\t')
                gid = line_split[0]
                gtdb_taxa = [x.strip() for x in line_split[1].split(';')]
                gtdb_taxonomy_from_file[gid] = gtdb_taxa

        Metadata = namedtuple('Metadata', ['ncbi_ref_genomes',
                                            'ncbi_rep_genomes',
                                            'ncbi_complete_genomes',
                                            'gtdb_representative',
                                            'genome_quality',
                                            'genome_comp',
                                            'genome_cont',
                                            'gtdb_taxonomy',
                                            'ncbi_taxonomy',
                                            'ncbi_org_name',
                                            'scaffold_count',
                                            'contig_count',
                                            'n50_scaffolds',
                                            'ambiguous_bases',
                                            'total_gap_length'])

        metadata = {}

        csv_reader = csv.reader(open(metadata_file))
        bHeader = True
        for row in csv_reader:
            if bHeader:
                genome_index = row.index('accession')
                org_index = row.index('organism_name')
                comp_index = row.index('checkm_completeness')
                cont_index = row.index('checkm_contamination')
                gtdb_taxonomy_index = row.index('gtdb_taxonomy')
                ncbi_taxonomy_index = row.index('ncbi_taxonomy')
                ncbi_org_name_index = row.index('ncbi_organism_name')
                rep_index = row.index('ncbi_refseq_category')
                assembly_level_index = row.index('ncbi_assembly_level')
                gtdb_rep_index = row.index('gtdb_representative')
                
                scaffold_count_index = row.index('scaffold_count')
                contig_count_index = row.index('contig_count')
                n50_scaffolds_index = row.index('n50_scaffolds')
                ambiguous_bases_index = row.index('ambiguous_bases')
                total_gap_length_index = row.index('total_gap_length')
                
                bHeader = False
            else:
                genome_id = row[genome_index]
                organism_name = row[org_index]
                if genome_id.startswith('U_') and '(UBA' not in organism_name:
                    continue
                
                comp = float(row[comp_index])
                cont = float(row[cont_index])
                quality = comp-5*cont
                
                scaffold_count = int(row[scaffold_count_index])
                contig_count = int(row[contig_count_index])
                n50_scaffolds = int(row[n50_scaffolds_index])
                ambiguous_bases = int(row[ambiguous_bases_index])
                total_gap_length = int(row[total_gap_length_index])
                                
                if gtdb_taxonomy_from_file:
                    if genome_id not in gtdb_taxonomy_from_file:
                        continue
                    gtdb_taxonomy = gtdb_taxonomy_from_file[genome_id]
                else:
                    gtdb_taxonomy = [x.strip() for x in row[gtdb_taxonomy_index].split(';')]
                ncbi_taxonomy = [x.strip() for x in row[ncbi_taxonomy_index].split(';')]
                ncbi_org_name = row[ncbi_org_name_index]
                ncbi_ref = row[rep_index]
                assembly_level = row[assembly_level_index]

                ncbi_ref_genomes = False
                ncbi_rep_genomes = False
                if ncbi_ref.strip() == 'reference genome':
                    ncbi_ref_genomes = True
                elif ncbi_ref.strip() == 'representative genome':
                    ncbi_rep_genomes = True
                    
                gtdb_representative = False
                if row[gtdb_rep_index] == 't':
                    gtdb_representative = True
                    
                ncbi_complete_genomes = False
                if assembly_level.strip() == 'Complete Genome':
                    ncbi_complete_genomes = True
                    
                metadata[genome_id] = Metadata(ncbi_ref_genomes=ncbi_ref_genomes,
                                                ncbi_rep_genomes=ncbi_rep_genomes,
                                                ncbi_complete_genomes=ncbi_complete_genomes,
                                                gtdb_representative=gtdb_representative,
                                                genome_quality=quality,
                                                genome_comp=comp,
                                                genome_cont=cont,
                                                gtdb_taxonomy=gtdb_taxonomy,
                                                ncbi_taxonomy=ncbi_taxonomy,
                                                ncbi_org_name=ncbi_org_name,
                                                scaffold_count=scaffold_count,
                                                contig_count=contig_count,
                                                n50_scaffolds=n50_scaffolds,
                                                ambiguous_bases=ambiguous_bases,
                                                total_gap_length=total_gap_length)   
                                                              
        return metadata
   
    def _genomes_in_rank(self, domain, metadata):
        """Get genomes in each taxa at each rank."""
        
        genomes_in_rank = defaultdict(lambda: defaultdict(set))
        for genome_id, m in metadata.items():
            taxa = m.gtdb_taxonomy
            
            if domain == 'archaea' and taxa[0] != 'd__Archaea':
                continue
            
            if domain == 'bacteria' and taxa[0] != 'd__Bacteria':
                continue
                
            for rank_index, taxon in enumerate(taxa):
                genomes_in_rank[rank_index][taxon].add(genome_id)
        
        return genomes_in_rank
        
    def _filter_genomes(self, 
                            genome_ids,
                            metadata,
                            quality_threshold,
                            canonical_count, min_canonical,
                            count_rps16, min_rps16,
                            count_rps23, min_rps23):
        """Filter genomes based on selection criteria."""
        
        quality_genome_ids = {}
        filtered_desc = {}
        for genome_id in genome_ids:
            qual = metadata[genome_id].genome_quality
            if qual < quality_threshold:
                filtered_desc[genome_id] = 'Failed quality test: %.1f' % qual
                continue
                
            if metadata[genome_id].scaffold_count > self.max_scaffolds:
                filtered_desc[genome_id] = 'Failed scaffold count test: %d' % metadata[genome_id].scaffold_count
                continue
                
            if metadata[genome_id].contig_count > self.max_contigs:
                filtered_desc[genome_id] = 'Failed contig count test: %d' % metadata[genome_id].contig_count
                continue
                
            if metadata[genome_id].n50_scaffolds < self.min_n50_scaffolds:
                filtered_desc[genome_id] = 'Failed scaffold N50 test: %d' % metadata[genome_id].n50_scaffolds
                continue
                
            if metadata[genome_id].ambiguous_bases > self.max_ambiguous:
                filtered_desc[genome_id] = 'Failed ambiguous bases test: %d' % metadata[genome_id].ambiguous_bases
                continue
                
            if metadata[genome_id].total_gap_length > self.max_total_gap_length:
                filtered_desc[genome_id] = 'Failed total gap length test: %d' % metadata[genome_id].total_gap_length
                continue
                
            if canonical_count.get(genome_id, -1) < min_canonical:
                filtered_desc[genome_id] = 'Failed ar122/bac120 test: %d' % canonical_count.get(genome_id, -1)
                continue
                  
            if count_rps16.get(genome_id, -1) < min_rps16:
                filtered_desc[genome_id] = 'Failed rps16 test: %d' % count_rps16.get(genome_id, -1)
                continue
                
            if count_rps23.get(genome_id, -1) < min_rps23:
                filtered_desc[genome_id] = 'Failed rps23 test: %d' % count_rps23.get(genome_id, -1)
                continue
                
            quality_genome_ids[genome_id] = qual
            
        return quality_genome_ids, filtered_desc
        
    def _ssu_test(self, 
                    genome_ids, metadata,
                    ssu_length, min_ssu_len):
        """Filter genomes based on selection criteria."""
        
        passed_ssu_test = set()
        for genome_id in genome_ids:
            if ssu_length.get(genome_id, -1) < min_ssu_len:
                continue

            passed_ssu_test.add(genome_id)
            
        return passed_ssu_test
        
    def _qualitative_metadata(self, genome_ids, metadata):
        """Get qualitative annotation information about genomes."""
        
        gtdb_rep = set()
        ncbi_ref = set()
        ncbi_rep = set()
        ncbi_complete = set()
        
        for genome_id in genome_ids:
            m = metadata[genome_id]
            
            if m.gtdb_representative:
                gtdb_rep.add(genome_id)
                
            if m.ncbi_ref_genomes:
                ncbi_ref.add(genome_id)
                
            if m.ncbi_rep_genomes:
                ncbi_rep.add(genome_id)
                
            if m.ncbi_complete_genomes:
                ncbi_complete.add(genome_id)
                
        return gtdb_rep, ncbi_ref, ncbi_rep, ncbi_complete
        
    def _sample_genomes(self,
                        genomes_per_taxon,
                        domain, 
                        rank_to_sample, 
                        genomes_in_rank, 
                        metadata,
                        quality_threshold,
                        ssu_length, min_ssu_len,
                        canonical_count, min_canonical,
                        count_rps16, min_rps16,
                        count_rps23, min_rps23,
                        output_dir):
        """Sample genomes from taxa."""
        
        rank_label = Taxonomy.rank_labels[rank_to_sample]
        print(rank_label)
        
        # sample genomes
        output_file = os.path.join(output_dir, 'selected_genomes.%s.%s.tsv' % (domain, rank_label))
        fout = open(output_file, 'w')
        fout.write('#Genome ID\t%s\tGTDB taxonomy\tNCBI taxonomy\tNCBI Organism Name\tCompleteness\tContamination\tQuality\tDescription\n' % rank_label.capitalize())
        
        missing_file = os.path.join(output_dir, 'unrepresented_taxa.%s.%s.tsv' % (domain, rank_label))
        fout_missing = open(missing_file, 'w')
        
        no_taxon_rep = set()
        ssu_status = {}
        for taxon, genome_ids in genomes_in_rank[rank_to_sample].items():
            # filter genomes based on general selection criteria
            quality_genome_ids, filtered_desc = self._filter_genomes(genome_ids,
                                                        metadata,
                                                        quality_threshold,
                                                        canonical_count, min_canonical,
                                                        count_rps16, min_rps16,
                                                        count_rps23, min_rps23)
                                                        
            # filter genomes based on 16S rRNA
            pass_ssu_test = self._ssu_test(quality_genome_ids.keys(),
                                            metadata, ssu_length, min_ssu_len)
                                            
            if pass_ssu_test.intersection(quality_genome_ids):
                # restrict genomes to those that passed filtering and have a SSU sequences
                failed_ssu_test = set(quality_genome_ids).difference(pass_ssu_test)
                for genome_id in failed_ssu_test:
                    quality_genome_ids.pop(genome_id)
                                  
                ssu_status[taxon] = True
            else:
                # based selection of genomes passing general filtering criteria and accept
                # that this taxon doesn't have a quality genome with a SSU sequence
                if len(quality_genome_ids): # otherwise the failure is due to other filtering
                    ssu_status[taxon] = False 
                                                        
            # get qualitative information about genomes
            gtdb_rep, ncbi_ref, ncbi_rep, ncbi_complete = self._qualitative_metadata(quality_genome_ids,
                                                                                        metadata)
            
            gtdb_rep_ids = gtdb_rep.intersection(quality_genome_ids)
            
            comp_ref_sp = gtdb_rep_ids.intersection(ncbi_complete).intersection(ncbi_ref)
            ref_sp = gtdb_rep_ids.intersection(ncbi_ref).difference(comp_ref_sp)
            
            comp_rep_sp = gtdb_rep_ids.intersection(ncbi_complete).intersection(ncbi_rep)
            rep_sp = gtdb_rep_ids.intersection(ncbi_rep).difference(comp_rep_sp)
            
            sampled_genomes = []
            if comp_ref_sp:
                quality = {genome_id:quality_genome_ids[genome_id] for genome_id in comp_ref_sp}
                quality_sorted = sorted(list(quality.items()), key=operator.itemgetter(1,0), reverse=True)
                sampled_genomes.extend([x[0] for x in quality_sorted[0:genomes_per_taxon - len(sampled_genomes)]])

            if ref_sp and len(sampled_genomes) != genomes_per_taxon:
                quality = {genome_id:quality_genome_ids[genome_id] for genome_id in ref_sp}
                quality_sorted = sorted(list(quality.items()), key=operator.itemgetter(1,0), reverse=True)
                sampled_genomes.extend([x[0] for x in quality_sorted[0:genomes_per_taxon - len(sampled_genomes)]])
            
            if comp_rep_sp and len(sampled_genomes) != genomes_per_taxon:
                quality = {genome_id:quality_genome_ids[genome_id] for genome_id in comp_rep_sp}
                quality_sorted = sorted(list(quality.items()), key=operator.itemgetter(1,0), reverse=True)
                sampled_genomes.extend([x[0] for x in quality_sorted[0:genomes_per_taxon - len(sampled_genomes)]])
                
            if rep_sp and len(sampled_genomes) != genomes_per_taxon:
                quality = {genome_id:quality_genome_ids[genome_id] for genome_id in rep_sp}
                quality_sorted = sorted(list(quality.items()), key=operator.itemgetter(1,0), reverse=True)
                sampled_genomes.extend([x[0] for x in quality_sorted[0:genomes_per_taxon - len(sampled_genomes)]])
                
            if gtdb_rep_ids and len(sampled_genomes) != genomes_per_taxon:
                quality = {genome_id:quality_genome_ids[genome_id] for genome_id in gtdb_rep_ids}
                quality_sorted = sorted(list(quality.items()), key=operator.itemgetter(1,0), reverse=True)
                sampled_genomes.extend([x[0] for x in quality_sorted[0:genomes_per_taxon - len(sampled_genomes)]])
                
            if len(sampled_genomes) != genomes_per_taxon:
                remaining_genomes = set(quality_genome_ids.keys()).difference(sampled_genomes)

                if remaining_genomes:
                    quality = {genome_id:quality_genome_ids[genome_id] for genome_id in remaining_genomes}
                    quality_sorted = sorted(list(quality.items()), key=operator.itemgetter(1,0), reverse=True)
                    sampled_genomes.extend([x[0] for x in quality_sorted[0:genomes_per_taxon - len(sampled_genomes)]])
                
            for genome_id in sampled_genomes:
                desc = 'Selected based on genome quality'
                if genome_id in comp_ref_sp:
                    desc = 'Complete NCBI reference genome'
                elif genome_id in ref_sp:
                    desc = 'NCBI reference genome'
                elif genome_id in comp_rep_sp:
                    desc = 'Complete NCBI representative genome'
                elif genome_id in rep_sp:
                    desc = 'NCBI representative genome'
                elif genome_id in gtdb_rep_ids:
                    desc = 'GTDB representative genome'
                
                fout.write('%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%s\n' % (genome_id, 
                                                                            taxon, 
                                                                            ';'.join(metadata[genome_id].gtdb_taxonomy), 
                                                                            ';'.join(metadata[genome_id].ncbi_taxonomy), 
                                                                            metadata[genome_id].ncbi_org_name,
                                                                            metadata[genome_id].genome_comp, 
                                                                            metadata[genome_id].genome_cont,
                                                                            metadata[genome_id].genome_quality,                                                               
                                                                            desc))
                                                                            
            if not sampled_genomes:
                fout_missing.write('%s\n' % taxon)
                for genome_id in genome_ids:
                    fout_missing.write('\t%s\t%s\n' % (genome_id, filtered_desc[genome_id]))
                no_taxon_rep.add(taxon)
                
        fout.close()
        fout_missing.close()
        
        if no_taxon_rep:
            print('%s: no reps for %d taxa' % (rank_label, len(no_taxon_rep)))
            
        print('%s: %d of %d taxa are missing a suitable 16S rRNA gene\n' % (rank_label, list(ssu_status.values()).count(False), len(ssu_status)))
        
    def run(self, 
            metadata_file,
            ar122_count_file,
            bac120_count_file,
            rps23_count_file, 
            rps16_ar_count_file, 
            rps16_bac_count_file,
            ssu_fasta_file,
            output_dir,
            genomes_per_taxon,
            taxonomy_file,
            quality_threshold,
            min_ssu_len,
            min_ar122,
            min_bac120,
            min_rps16,
            min_rps23):
        """Sample taxa."""
        
        # create output directory
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # record parameters
        param_file = os.path.join(output_dir, 'dereplicate_taxa.log')
        fout = open(param_file, 'w')
        fout.write(ntpath.basename(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]) + '\n')
        fout.close()
        
        # get genome metadata
        print('Reading metadata.')
        metadata = self._read_metadata(metadata_file, taxonomy_file)
        
        # get count for ribosomal marker sets
        print('Reading count of each marker set.')
        count_ar122 = self._read_marker_count_file(ar122_count_file)
        count_bac120 = self._read_marker_count_file(bac120_count_file)
        count_rps23 = self._read_marker_count_file(rps23_count_file)
        count_rps16_ar = self._read_marker_count_file(rps16_ar_count_file)
        count_rps16_bac = self._read_marker_count_file(rps16_bac_count_file)
        
        # get length of SSU sequences
        print('Reading SSU length.')
        ssu_length = self._read_ssu_file(ssu_fasta_file)

        # create combined rps16 set
        count_rps16 = {}
        for genome_id, m in metadata.items():
            gtdb_domain = m.gtdb_taxonomy[0]
            if gtdb_domain == 'd__Bacteria':
                count_rps16[genome_id] = count_rps16_bac.get(genome_id, -1)
            elif gtdb_domain == 'd__Archaea':
                count_rps16[genome_id] = count_rps16_ar.get(genome_id, -1)
            else:
                count_rps16[genome_id] = max(count_rps16_bac.get(genome_id, -1), 
                                                count_rps16_ar.get(genome_id, -1))
            
        # sample 1 genome per taxonomic group at the ranks of phylum to genus
        print('Selecting representatives.')
        genomes_per_rank = 1
        for domain in ['archaea', 'bacteria']:
            print(domain)
            
            # get canonical marker set information
            if domain == 'archaea':
                canonical_count = count_ar122
                min_canonical = min_ar122
            else:
                canonical_count = count_bac120
                min_canonical = min_bac120
            
            # get genomes for each taxon
            genomes_in_rank = self._genomes_in_rank(domain, metadata)
            
            # identify genera without any named species
            genera_without_named_species = 0
            for genus, genome_ids in genomes_in_rank[5].items():
                had_named_species = False
                for genome_id in genome_ids:
                    if metadata[genome_id].gtdb_taxonomy[6] != 's__':
                        had_named_species = True
                        break
                
                if not had_named_species:
                    genera_without_named_species += 1
                    genus_taxon = genus.replace('g__', '')
                    genomes_in_rank[6]['s__' + genus_taxon + '_unclassified'] = genome_ids
                    
            print('Identified %d genera without a named species.' % genera_without_named_species)
            print('For species dereplication, a genome will be selected from the genus.')
            
            # dereplicate from each taxon at each rank
            for rank_to_sample in range(1, 7):
                self._sample_genomes(genomes_per_taxon,
                                        domain, 
                                        rank_to_sample, 
                                        genomes_in_rank, 
                                        metadata,
                                        quality_threshold,
                                        ssu_length, min_ssu_len,
                                        canonical_count, min_canonical,
                                        count_rps16, min_rps16,
                                        count_rps23, min_rps23,
                                        output_dir)

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('metadata_file', help='metadata for GTDB genomes')
    parser.add_argument('ar122_count_file', help='count information for ar122 markers (determined with genedb)')
    parser.add_argument('bac120_count_file', help='count information for bac120 markers (determined with genedb)')
    parser.add_argument('rps23_count_file', help='count information for rps23 markers (determined with genedb)')
    parser.add_argument('rps16_ar_count_file', help='count information for archaeal rps16 markers (determined with genedb)')
    parser.add_argument('rps16_bac_count_file', help='count information for bacterial rps16 markers (determined with genedb)')
    parser.add_argument('ssu_fasta_file', help='FASTA file with SSU sequences from GTDB')
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('--genomes_per_taxon', help='genomes to sample from each taxon', type=int, default=1)

    parser.add_argument('--taxonomy', help='taxonomy to use instead of GTDB taxonomy specified in the metadata file')
    parser.add_argument('--quality_threshold', help='minimum quality to select genome', type=float, default=50.0)
    parser.add_argument('--min_ssu_len', help='minimum length of 16S rRNA gene', type=int, default=900)
    parser.add_argument('--min_ar122', help='minimum number of single-copy genes in ar122 marker set', type=int, default=61) 
    parser.add_argument('--min_bac120', help='minimum number of single-copy genes in bac120 marker set', type=int, default=60) 
    parser.add_argument('--min_rps16', help='minimum number of single-copy genes in rps16 marker set (18 HMMs)', type=int, default=9)
    parser.add_argument('--min_rps23', help='minimum number of single-copy genes in rps23 marker set (27 HMMs)', type=int, default=14)
    
    parser.add_argument('--min_n50', help='minimum scaffold N50 to consider genome', type=int, default=10000)
    parser.add_argument('--max_contig_count', help='maximum contig count to consider genome', type=int, default=1000)
    parser.add_argument('--max_scaffold_count', help='maximum scaffolds count to consider genome', type=int, default=500)
    
    args = parser.parse_args()

    try:
        p = DereplicateTaxa(args.min_n50, args.max_contig_count, args.max_scaffold_count)
        p.run(args.metadata_file,
                args.ar122_count_file,
                args.bac120_count_file,
                args.rps23_count_file, 
                args.rps16_ar_count_file, 
                args.rps16_bac_count_file,
                args.ssu_fasta_file,
                args.output_dir,
                args.genomes_per_taxon,
                args.taxonomy,
                args.quality_threshold,
                args.min_ssu_len,
                args.min_ar122,
                args.min_bac120,
                args.min_rps16,
                args.min_rps23)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
