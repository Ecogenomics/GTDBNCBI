# import prodigal
import multiprocessing
import math
import time
import random

from itertools import islice
from gtdblite.Exceptions import GenomeDatabaseError
from gtdblite import MarkerCalculation
from itertools import islice


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
