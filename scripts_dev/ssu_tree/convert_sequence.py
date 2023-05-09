from Bio import SeqIO
input_file = open("ssu_all_r214.fna")
seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(input_file, "fasta")}

my_dict = {k.split('~')[0]:v for k,v in seq_dict.items()}

list_reps = []
with open('bac120_reps.lst') as gtb:
    for line in gtb:
        list_reps.append(line.strip())

fasta_aligned = open('reps_ssu_non_aligned_gt700.fna','w')
for rep in list_reps:
    if rep in my_dict and len(my_dict.get(rep))>=700:
        fasta_aligned.write(f'>{rep}\n{my_dict.get(rep)}\n')
fasta_aligned.close()
