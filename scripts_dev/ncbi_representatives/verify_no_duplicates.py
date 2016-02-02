# Verify that no genome is clustered multiple time.

import sys

s = set()
for line in open('gtdb_clusters_99.8.tsv'):
    line_split = line.strip().split('\t')
    
    rep_id = line_split[0]
    if rep_id in s:
        print 'Duplicate: ', s
    s.add(rep_id)
    
    if len(line_split) == 4:
        genome_ids = line_split[3].split(',')
        for genome_id in genome_ids:
            if genome_id in s:
                print 'Duplicate: ', s
                sys.exit()
                
            s.add(genome_id)
    
print len(s)
