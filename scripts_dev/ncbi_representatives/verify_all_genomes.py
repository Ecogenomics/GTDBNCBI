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

            genome_quality[genome_id] = [comp, cont, comp - 4*cont]

    return genome_quality

reps = {}
for line in open('representatives.tsv'):
    line_split = line.strip().split('\t')
    rep = line_split[0]
    
    reps[rep] = set([rep])
    
    if len(line_split) == 4:
        genomes = set(line_split[3].split(','))
        reps[rep].update(genomes)
    
print 'Reps: %d' % len(reps)
print 'Genomes: %d' % sum([len(x) for x in reps.values()])

genome_ids = set()
for line in open('genomes_to_retain.tsv'):
    line_split = line.split('\t')
    genome_id = line_split[0]
    
    if genome_id in reps:
        genome_ids.update(reps[genome_id])
    else:
        genome_ids.add(genome_id)
        
print 'Total genomes retained or represented: %d' % len(genome_ids)

genome_quality = read_gtdb_genome_quality('../../gtdb_metadata.csv')
print 'Total genomes in database: %d' % len(genome_quality)

for g in genome_quality:
    if g not in genome_ids:
        comp, cont, qual = genome_quality[g]
        if comp >= 10.0 and cont <= 10.0 and qual > 0:
            print g, genome_quality[g]