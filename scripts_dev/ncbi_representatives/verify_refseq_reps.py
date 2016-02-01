# Verify that initial set of RefSeq representatives
# are retained in final set of representatives.

refseq_reps = set()
for line in open('../../ncbi_refseq_representatives.tsv'):
    if line[0] == '#':
        continue
        
    line_split = line.split('\t')
    refseq_reps.add(line_split[0])
    
for line in open('../gtdb_clusters_99.tsv'):
    line_split = line.split('\t')
    
    rep = line_split[0]
    if rep in refseq_reps:
        refseq_reps.remove(rep)
        
print len(refseq_reps)