s = set()
for line in open('genomes_to_retain.tsv'):
    line_split = line.split('\t')
    genome_id = line_split[0]
    
    if genome_id in s:
        print 'Duplicate: %s' % genome_id
        
    s.add(genome_id)
    
print len(s)
