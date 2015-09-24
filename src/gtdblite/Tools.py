import prodigal
import tempfile
import os
import shutil
import re
import subprocess
import multiprocessing
import math
import time

from gtdblite import GenomeDatabase
from gtdblite import profiles
from gtdblite.Exceptions import GenomeDatabaseError

##################################################
############FILE UTILITIES########################
##################################################


def populate_required_headers(checkm_fh):
    required_headers = {
	"Bin Id" : None,
	"Completeness" : None,
	"Contamination" : None
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



##################################################
############MISC UTILITIES########################
##################################################
def splitchunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def generateTempTableName(self):
    rng = random.SystemRandom()
    suffix = ''
    for i in range(0,10):
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
    
