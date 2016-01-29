import csv

def read_gtdb_taxonomy(metadata_file):
    """Parse GTDB taxonomy from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes.

    Return
    ------
    dict : d[genome_id] -> taxonomy list
    """

    taxonomy = {}

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('genome')
            taxonomy_index = headers.index('gtdb_taxonomy')
            bHeader = False
        else:
            genome_id = row[genome_index]
            taxonomy[genome_id] = row[taxonomy_index].split(';')

    return taxonomy


def read_gtdb_ncbi_taxonomy(metadata_file):
    """Parse NCBI taxonomy from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes.

    Return
    ------
    dict : d[genome_id] -> taxonomy list
    """

    taxonomy = {}

    csv_reader = csv.reader(open(metadata_file, 'rt'))
    bHeader = True
    for row in csv_reader:
        if bHeader:
            headers = row
            genome_index = headers.index('genome')
            taxonomy_index = headers.index('ncbi_taxonomy')
            bHeader = False
        else:
            genome_id = row[genome_index]
            taxonomy[genome_id] = row[taxonomy_index].split(';')

    return taxonomy
    
taxonomy = read_gtdb_taxonomy('../../gtdb_metadata.csv')
#taxonomy = read_gtdb_ncbi_taxonomy('../../gtdb_metadata.csv')

for line in open('representatives.tsv'):
    line_split = line.strip().split('\t')
    
    rep_id = line_split[0]
    rep_taxonomy = taxonomy[rep_id]
    
    if len(rep_taxonomy) == 7 and len(line_split) == 4:
        for genome_id in line_split[3].split(','):
            genome_taxonomy = taxonomy[genome_id]
            if len(genome_taxonomy) == 7:
                rep_sp = rep_taxonomy[6]
                genome_sp = genome_taxonomy[6]
                
                if rep_sp != genome_sp and genome_sp != 's__' and rep_sp != 's__':
                    print 'Inter-species clustering:', rep_taxonomy[5], rep_sp, genome_sp
    