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
    
for line in open('ani_99.8.ani.tsv'):
    line_split = line.strip().split('\t')
    rep_id = line_split[0]
    genome_id = line_split[1]
    ani = float(line_split[2])
    
    if ani < 96.5 and ani > 0.0 and genome_quality[genome_id][2] > 80:
        print line_split, genome_quality[rep_id], genome_quality[genome_id]