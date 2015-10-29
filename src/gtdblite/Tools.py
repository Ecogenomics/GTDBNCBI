# import prodigal
import hashlib
import tempfile
import os
import shutil
import re
import subprocess
import multiprocessing
import math
import time
import random

from itertools import islice
from gtdblite.Exceptions import GenomeDatabaseError
from gtdblite import MarkerCalculation
from gtdblite import Config

##################################################
##################################################
############FILE UTILITIES########################
##################################################


def populate_required_headers(checkm_fh):
    required_headers = {
        "Bin Id": None,
        "Completeness": None,
        "Contamination": None,
        "Marker lineage": None,
        "# genomes": None,
        "# markers": None,
        "# marker sets": None,
        "Strain heterogeneity": None
    }

    # Check the CheckM headers are consistent
    split_headers = checkm_fh.readline().rstrip().split("\t")
    for pos in range(0, len(split_headers)):
        header = split_headers[pos]
        if header not in required_headers:
            continue

        if required_headers[header] is not None:
            raise GenomeDatabaseError(
                "Seen %s header twice in the CheckM file. Check that the CheckM file is correct: %s." % (header, checkm_fh.name))

        required_headers[header] = pos

    for header, col in required_headers.items():
        if (header is "Completeness" or header is "Contamination") and col is None:
            raise GenomeDatabaseError(
                "Unable to find %s header in the CheckM file. Check that the CheckM file is correct: %s." % (header, checkm_fh.name))

    return required_headers


#################################################
#################################################
#################################################
def runMultiProdigal(nprocs=None, nogene_list=None):
    def worker(nogene_subdict, out_dict):
        """This worker function is invoked in a process."""

        for key, value in nogene_subdict.iteritems():
            if value.get("aa_gene_file") is None:
                rtn_files = MarkerCalculation.RunProdigalOnGenomeFasta(
                    value.get("fasta_path"))
                aa_gene_file, nt_gene_file, gff_file, translation_table_file = rtn_files
                value["aa_gene_file"] = aa_gene_file
                value["nt_gene_file"] = nt_gene_file
                value["gff_file"] = gff_file
                value["translation_table_file"] = translation_table_file

            out_dict[key] = value

        return True

    # split bins to process over specified number of processors
    manager = multiprocessing.Manager()
    out_dict = manager.dict()

    procs = []
    for item in splitchunks(nogene_list, nprocs):
        p = multiprocessing.Process(target=worker, args=(item, out_dict))
        procs.append(p)
        p.start()

    # Collect all results into a single result dict. We know how many dicts
    # with results to expect.
    while len(out_dict) < len(nogene_list):
        time.sleep(1)

    # Wait for all worker processes to finish
    for p in procs:
        p.join()

    return out_dict


##################################################
############MISC UTILITIES########################
##################################################
def splitchunks(d, n):
    chunksize = int(math.ceil(len(d) / float(n)))
    it = iter(d)
    for _ in xrange(0, len(d), chunksize):
        yield {k: d[k] for k in islice(it, chunksize)}


def generateTempTableName():
    rng = random.SystemRandom()
    suffix = ''
    for _ in range(0, 10):
        suffix += rng.choice(
            'abcefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    return "TEMP" + suffix + str(int(time.time()))


def fastaPathGenerator(path=None, prefix=None):

    genomeUserDir = None
    if Config.GTDB_GENOME_USR_DIR:
        genomeUserDir = Config.GTDB_GENOME_USR_DIR

    genomeNCBIDir = None
    if Config.GTDB_GENOME_NCBI_DIR:
        genomeNCBIDir = Config.GTDB_GENOME_NCBI_DIR

    if prefix is 'U':
        return os.path.join(genomeUserDir, path)
    elif prefix is 'NCBI':
        return os.path.join(genomeNCBIDir, path)
