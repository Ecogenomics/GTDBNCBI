# import prodigal
import math
import time
import random
import os

from itertools import islice

import Config


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

    genomeGBKDir = None
    if Config.GTDB_GENOME_GBK_DIR:
        genomeGBKDir = Config.GTDB_GENOME_GBK_DIR

    genomeRSQDir = None
    if Config.GTDB_GENOME_RSQ_DIR:
        genomeRSQDir = Config.GTDB_GENOME_RSQ_DIR

    if prefix == 'U':
        return os.path.join(genomeUserDir, path)
    elif prefix == "GB":
        return os.path.join(genomeGBKDir, path)
    elif prefix == "RS":
        return os.path.join(genomeRSQDir, path)
    else:
        print "prefix {0} is not existing".format(prefix)


def confirm(msg):
    raw = raw_input(msg + " (y/N): ")
    if raw.upper() == "Y":
        return True
    return False
