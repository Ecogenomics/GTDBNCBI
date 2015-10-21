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
from itertools import islice

##################################################
##################################################
############FILE UTILITIES########################
##################################################


def populate_required_headers(checkm_fh):
    required_headers = {
    "Bin Id": None,
    "Completeness": None,
    "Contamination": None
    }

    # Check the CheckM headers are consistent
    split_headers = checkm_fh.readline().rstrip().split("\t")
    for pos in range(0, len(split_headers)):
        header = split_headers[pos]
        if header not in required_headers:
            continue
        if required_headers[header] is not None:
            raise GenomeDatabaseError("Seen %s header twice in the checkM file. Check the checkM file is correct: %s." % (header, checkM_file))
        required_headers[header] = pos
        for header, col in required_headers.items():
            if col is None:
                raise GenomeDatabaseError("Unable to find %s header in the checkM file. Check the checkM file is correct: %s." % (header, checkM_file))
    return required_headers


#################################################
#################################################
#################################################
def runMultiProdigal(nprocs=None, nogene_list=None):
    def worker(nogene_subdict, out_dict):
        """ The worker function, invoked in a process. 'genome_dict' is a
            dictionnary of genes to align. The results are placed in
            a dictionary that's pushed to a queue.
        """

        for key, value in nogene_subdict.iteritems():
            if(value.get("gene_path") is None):
                prodigal_tmp_dir = MarkerCalculation.RunProdigalOnGenomeFasta(value.get("fasta_path"))
                value["gene_path"] = os.path.join(prodigal_tmp_dir, "genes.faa")
            out_dict[key] = value
        return True

    # Each process will get 'chunksize' nums and a queue to put his out
    # dict into
    manager = multiprocessing.Manager()
    out_dict = manager.dict()
#    chunksize = int(math.ceil(len(genome_dict) / float(nprocs)))
    procs = []

    for item in splitchunks(nogene_list, nprocs):
        p = multiprocessing.Process(
                target=worker,
                args=(item, out_dict))
        procs.append(p)
        p.start()

    # Collect all results into a single result dict. We know how many dicts
    # with results to expect.
    while len(out_dict) < len(nogene_list):
        time.sleep(1)

    # Wait for all worker processes to finish
    for p in procs:
        p.join()

    print out_dict

    return out_dict


##################################################
############MISC UTILITIES########################
##################################################
def splitchunks(d, n):
    chunksize = int(math.ceil(len(genome_dict) / float(nprocs)))
    it = iter(d)
    for _ in xrange(0, len(d), chunksize):
        yield {k: d[k] for k in islice(it, chunksize)}


def generateTempTableName():
    rng = random.SystemRandom()
    suffix = ''
    for _ in range(0, 10):
        suffix += rng.choice('abcefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    return "TEMP" + suffix + str(int(time.time()))


def sha256Calculator(file_path):
    try:
        filereader = open(file_path, "rb")
    except:
        raise GenomeDatabaseError("Cannot open Fasta file: " + file_path)

    m = hashlib.sha256()
    for line in filereader:
        m.update(line)
    sha256_checksum = m.hexdigest()
    filereader.close()
    return sha256_checksum
