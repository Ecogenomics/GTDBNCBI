# import prodigal
import math
import time
import random
import os

from itertools import islice
from gtdblite.Exceptions import GenomeDatabaseError
from gtdblite import MarkerCalculation
from gtdblite import Config
from mhlib import PATH

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

    if prefix == 'U':
        return os.path.join(genomeUserDir, path)
    elif prefix == "NCBI":
        return os.path.join(genomeNCBIDir, path)
    else:
        print "prefix {0} is not existing".format(prefix)
