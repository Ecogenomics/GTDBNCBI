# Verify that all representative genomes have
# the expected genome quality (completeness - contamination)

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

reps = set()
for line in open('../gtdb_representatives_99.txt'):
    line_split = line.strip().split('\t')
    rep = line_split[0]
    
    reps.add(rep)

genome_quality = read_gtdb_genome_quality('../../../gtdb_metadata.csv')

fout = open('rep_genome_quality.tsv', 'w')
for r in reps:
    comp, cont, qual = genome_quality[r]
    fout.write('%s\t%.2f\t%.2f\t%.2f\n' % (r, float(comp), float(cont), float(qual)))
fout.close()
    
