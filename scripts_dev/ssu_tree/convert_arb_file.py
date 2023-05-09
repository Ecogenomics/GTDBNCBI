import gzip
import os
import sys
import traceback


def canonical_gid(gid: str) -> str:
    """Get canonical form of NCBI genome accession.

    Example:
        G005435135 -> G005435135
        GCF_005435135.1 -> G005435135
        GCF_005435135.1_ASM543513v1_genomic -> G005435135
        RS_GCF_005435135.1 -> G005435135
        GB_GCA_005435135.1 -> G005435135
    """

    if gid.startswith('U'):
        return gid

    gid = gid.replace('RS_', '').replace('GB_', '')
    gid = gid.replace('GCA_', 'G').replace('GCF_', 'G')
    if '.' in gid:
        gid = gid[0:gid.find('.')]

    return gid


def read_fasta(fasta_file, keep_annotation=False):
    """Read sequences from fasta file.

    Parameters
    ----------
    fasta_file : str
        Name of fasta file to read.
    keep_annotation : boolean
        Determine is sequence id should contain annotation.

    Returns
    -------
    dict : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    """

    if not os.path.exists(fasta_file):
        raise Exception('Input file %s does not exist.' % fasta_file)

    if os.stat(fasta_file).st_size == 0:
        return {}

    try:

        if fasta_file.endswith('.gz'):
            file_f, file_mode = gzip.open, 'rt'
        else:
            file_f, file_mode = open, 'r'

        seqs = {}
        with file_f(fasta_file, file_mode) as f:

            for line in f.readlines():
                # skip blank lines
                if not line.strip():
                    continue

                if line[0] == '>':
                    if keep_annotation:
                        seq_id = line[1:-1]
                    else:
                        seq_id = line[1:].split(None, 1)[0]

                    seqs[seq_id] = []
                else:
                    seqs[seq_id].append(line.strip())

        for seq_id, seq in seqs.items():
            seqs[seq_id] = ''.join(seq).replace(' ', '')
    except:
        print(traceback.format_exc())
        print()
        print("[Error] Failed to process sequence file: " + fasta_file)
        sys.exit(1)

    return seqs

# we parse the fasta sequences

seqs= read_fasta('reps_ssu_aligned_4frames_NR99_gt700.fasta')

outf = open('converted_arb_ssu_gt700.txt','w')

with open('gtdb_arb_metadata.txt','r') as f:
    iteration =0
    keep = False
    for line in f.readlines():
        line = line.strip('\n')
        infos = line.split('=')

        if infos[0] == 'BEGIN':
            iteration +=1
            print(iteration)
        if infos[0] == 'db_name' and infos[1] in seqs:
            id_seq = infos[1]
            outf.write(f'BEGIN\n{infos[0]}={canonical_gid(infos[1])}\n')
            keep = True
        elif keep:
            if infos[0] == 'aligned_seq':
                outf.write(f'{infos[0]}={seqs.get(id_seq)}\n')
            elif infos[0] == 'END':
                outf.write(f'END\n\n')
                keep = False
            else:
                outf.write(f'{line}\n')