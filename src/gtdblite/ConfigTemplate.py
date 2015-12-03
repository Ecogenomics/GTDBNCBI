GTDB_HOST = "host"
GTDB_USERNAME = "dbuser"
GTDB_DB_NAME = "database_name"
GTDB_PASSWORD = "database_password"

# Root directory for genomes
GTDB_GENOME_ROOT = "/root/path/to/genomes/"

# Deprecated folder
GTDB_DPRCTD_ROOT = "/root/path/to/deprecated/genomes/"
USER_SUFFIX = "suffix/for/user/genomes/"
NCBI_SUFFIX = "suffix/for/ncbi/genomes/"
GBK_PREFIX = NCBI_SUFFIX + "genbank/"
RSQ_PREFIX = NCBI_SUFFIX + "refseq/"


GTDB_GENOME_USR_DIR = GTDB_GENOME_ROOT + USER_SUFFIX
GTDB_GENOME_GBK_DIR = GTDB_GENOME_ROOT + GBK_PREFIX
GTDB_GENOME_RSQ_DIR = GTDB_GENOME_ROOT + RSQ_PREFIX

GTDB_DPRCTD_USR_DIR = GTDB_DPRCTD_ROOT + USER_SUFFIX
GTDB_DPRCTD_GBK_DIR = GTDB_DPRCTD_ROOT + GBK_PREFIX
GTDB_DPRCTD_RSQ_DIR = GTDB_DPRCTD_ROOT + RSQ_PREFIX
