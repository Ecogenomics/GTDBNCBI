# Representatives are selected based on their
# source reporitory (RefSeq=>GenBank=>User),
# followed by genome quality. As such, there
# should be few, if any, cases were a representative
# User genome represents NCBI genomes.

import csv

def read_gtdb_genome_quality(metadata_file):
    genome_quality = {}

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('genome')
            comp_index = headers.index('checkm_completeness')
            cont_index = headers.index('checkm_contamination')
            bHeader = False
        else:
            genome_id = row[genome_index]
            comp = float(row[comp_index])
            cont = float(row[cont_index])

            genome_quality[genome_id] = [comp, cont, comp - cont]

    return genome_quality

genome_quality = read_gtdb_genome_quality('../../../gtdb_metadata.csv')


for line in open('gtdb_clusters_99.5.tsv'):
    line_split = line.strip().split('\t')
    ref_id = line_split[0]
    
    if ref_id.startswith('U_') and len(line_split) == 4:
        genome_ids = line_split[3].split(',')
        for genome_id in genome_ids:
            if not genome_id.startswith('U_'):
                print '%s\t%s\t%s\t%s' % (ref_id, str(genome_quality[ref_id]), genome_id, str(genome_quality[genome_id]))