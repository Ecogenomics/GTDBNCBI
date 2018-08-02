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

GTDB_GENOME_USR_DIR = GTDB_GENOME_ROOT + USER_SUFFIX
GTDB_GENOME_GBK_DIR = GTDB_GENOME_ROOT + NCBI_SUFFIX + 'genbank'
GTDB_GENOME_RSQ_DIR = GTDB_GENOME_ROOT + NCBI_SUFFIX + 'refseq'

GTDB_DPRCTD_USR_DIR = GTDB_DPRCTD_ROOT
GTDB_DPRCTD_GBK_DIR = GTDB_DPRCTD_ROOT
GTDB_DPRCTD_RSQ_DIR = GTDB_DPRCTD_ROOT

# Annotation folder
NCBI_ANNOTATION_DIR = 'prodigal'
USER_ANNOTATION_DIR = NCBI_ANNOTATION_DIR

# Block Insertion/Deletion of genomes during updates
# True : Stops Add/Delete
# False : Allows Add/Delete
DB_UPDATE = False

# List of available database and the server associated with each of them.
DB_SERVERS = {'deprecated_gtdb_version': 'hostname',
              'latest_gtdb_version': 'hostname'}
LATEST_DB = 'latest_gtdb_version'  # WHERE ADD AND DELETE CAN HAPPEN
