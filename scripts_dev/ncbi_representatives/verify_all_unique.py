genome_ids = set()
for line in open('representatives.tsv'):
    line_split = line.strip().split('\t')
    
    rep = line_split[0]
    if rep in genome_ids:
        print 'Duplicate rep genome: %s' % rep
        exit(-1)
        
    genome_ids.add(rep)
    
    if len(line_split) == 4:
        genomes = line_split[3].split(',')
        for g in genomes:
            if g in genome_ids:
                print 'Duplicate clustered genome: %s' % g
                exit(-1)
                
            genome_ids.add(g)
            
print 'Number of genomes: %d' % len(genome_ids)