import logging
import psycopg2

from gtdblite.GenomeDatabaseConnection import GenomeDatabaseConnection
from gtdblite.Exceptions import GenomeDatabaseError


class MetadataManager(object):

    def __init__(self):
        self.conn = GenomeDatabaseConnection()

    # Function: viewMetadata
    #
    # Lists all metadata field available in the database
    #
    # Returns:
    # Print lists
    def viewMetadata(self):
        '''
        # Function: viewMetadata
    #
    # Lists all metadata field available in the database
    #
    # Returns:
    # Print lists
        '''
        try:
            self.conn.MakePostgresConnection()
            cur = self.conn.cursor()
            cur.execute("SELECT * FROM view_list_meta_columns")
            print "\t".join(("Table", "Field", "Datatype", "Description"))
            for tup in cur.fetchall():
                info_list = list(tup)
                if info_list[2] == "double precision":
                    info_list[2] = "float"
                elif str(info_list[2]).startswith("timestamp"):
                    info_list[2] = "timestamp"
                print "\t".join(info_list)
            cur.close()
            self.conn.ClosePostgresConnection()
        except GenomeDatabaseError as e:
            raise e

    def exportMetadata(self, path):
        try:
            self.conn.MakePostgresConnection()
            cur = self.conn.cursor()
            query = "SELECT * from metadata_view"
            outputquery = 'copy ({0}) to stdout with csv header'.format(query)
            with open(path, 'w') as f:
                cur.copy_expert(outputquery, f)
            print "Export Successful"
            cur.close()
            self.conn.ClosePostgresConnection()
        except GenomeDatabaseError as e:
            raise e

    def importMetadata(self, table=None, field=None, typemeta=None, metafile=None):
        try:
            self.conn.MakePostgresConnection()
            cur = self.conn.cursor()
            data_list = []
            with open(metafile, 'r') as metaf:
                for line in metaf:
                    data_list.append(tuple(line.strip().split('\t')))
            data_zip = zip(*data_list)
            genome_id = list(data_zip[0])
            meta_value = list(data_zip[1])
            for n, i in enumerate(genome_id):
                new_i = i.split("_", 1)[1]
                genome_id[n] = new_i
            query = "SELECT upsert('{0}','{1}','{2}',%s,%s)".format(
                table, field, typemeta)
            cur.execute(query, (genome_id, meta_value))
            self.conn.commit()
            cur.close()
            self.conn.ClosePostgresConnection()
        except GenomeDatabaseError as e:
            raise e
        except psycopg2.Error as e:
            raise GenomeDatabaseError(e.pgerror)

    def createMetadata(self, metadatafile):
        try:
            self.conn.MakePostgresConnection()
            cur = self.conn.cursor()
            data_dict = {}
            with open(metadatafile, 'r') as metaf:
                for line in metaf:
                    array_line = line.strip().split('\t')
                    if not array_line[3].startswith("metadata_"):
                        raise GenomeDatabaseError(
                            "Only Metadata Tables can be modified")
                    data_dict[array_line[0]] = {
                        "table": array_line[3], "type": array_line[2], "desc": array_line[1]}

            query = "SELECT v.field,v.table from view_list_meta_columns as v"
            cur.execute(query)
            all_col_dict = dict(cur.fetchall())
            print all_col_dict
            for key, value in data_dict.iteritems():
                if key in all_col_dict:
                    if all_col_dict.get(key) == value.get("table"):
                        query_comment = "COMMENT ON COLUMN {0}.{1} IS '{2}'".format(
                            value.get("table"), key, value.get("desc"))
                        cur.execute(query_comment)
                    else:
                        logging.warning("Column {0} is already presents in the {1} table .".format(
                            key, all_col_dict.get(key)))
                else:
                    query_add_col = "ALTER TABLE {0} ADD COLUMN {1} {2}".format(
                        value.get("table"), key, value.get("type"))
                    cur.execute(query_add_col)
                    query_add_comment = "COMMENT ON COLUMN {0}.{1} IS '{2}'".format(
                        value.get("table"), key, value.get("desc"))
                    cur.execute(query_add_comment)

# ---------- PSQL are not refresh automatically so we need to drop the existing view and recreate it with a new Definition.
            cur.execute("SELECT refreshView()")
            self.conn.commit()
            cur.close()
            self.conn.ClosePostgresConnection
        except GenomeDatabaseError as e:
            raise e
