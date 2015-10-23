import os
import logging
import psycopg2
import sys


from psycopg2.extensions import AsIs

from gtdblite import ConfigMetadata
from gtdblite.GenomeDatabaseConnection import GenomeDatabaseConnection
from gtdblite.Exceptions import GenomeDatabaseError

from biolib.common import remove_extension


class MetadataManager(object):

    def __init__(self):
        self.conn = GenomeDatabaseConnection()
        self.conn.MakePostgresConnection()

    def viewMetadata(self):
        '''
         Function: viewMetadata
        Lists all metadata field available in the database

        Returns:
        Print lists of Metadata on Screen
        '''
        try:
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
            raise self.ReportError(e.message)

#     def exportMetadata(self, path):
#         try:
#             self.conn.MakePostgresConnection()
#             cur = self.conn.cursor()
#             query = "SELECT exportMeta('{0}')".format(path)
#             cur.execute(query)
#             print "Export Successful"
#             cur.close()
#             self.conn.ClosePostgresConnection
#         except GenomeDatabaseError as e:
#             raise self.ReportError(e.message)

    def exportMetadata(self, path):
        '''
        Function: exportMetadata
        Export all Metadata for all genomes to a csv file

        :param path: Path to the output file
        '''

        try:
            cur = self.conn.cursor()
            query = "SELECT * from metadata_view"
            outputquery = 'copy ({0}) to stdout with csv header'.format(query)
            with open(path, 'w') as f:
                cur.copy_expert(outputquery, f)
            print "Export Successful"
            cur.close()
            self.conn.ClosePostgresConnection()
        except GenomeDatabaseError as e:
            raise self.ReportError(e.message)

    def importMetadata(self, table=None, field=None, typemeta=None, metafile=None):
        '''
        Function importMetadata
        import one field of Metadata for a list of Genomes

        :param table: Table where the column is located
        :param field: Name of the Column
        :param typemeta: Data type of the column
        :param metafile: TSV file with the format (Genome_id \t Value)
        '''
        try:
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
            try:
                cur.execute(query, (genome_id, meta_value))
                self.conn.commit()
            except psycopg2.Error as e:
                print e.pgerror
            cur.close()
            self.conn.ClosePostgresConnection()
        except GenomeDatabaseError as e:
            raise self.ReportError(e.message)

    def createMetadata(self, metadatafile):
        '''
        Function createMetadata
        Create or Update metaddata columns in the database

        :param metadatafile: TSV file listing one new field per line
        Format of the TSV file is new_field \t description \t type \t table
        '''
        try:
            cur = self.conn.cursor()
            data_dict = {}
            with open(metadatafile, 'r') as metaf:
                for line in metaf:
                    array_line = line.strip().split('\t')
#                    if not array_line[3].startswith("metadata_"):
#                        raise GenomeDatabaseError(
#                            "Only Metadata Tables can be modified")
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
            self.conn.commit()
            cur.close()
            self.conn.ClosePostgresConnection
        except GenomeDatabaseError as e:
            raise self.ReportError(e.message)

    def calculateMetadata(self, genome_file, gff_file, output_dir):
        """Calculate metadata for new genome.

        Parameters
        ----------
        genome_file : str
            Name of fasta file containing nucleotide sequences.
        gff_file : str
            Name of generic feature file describing genes.
        output_dir : str
            Output directory.
        """

        os.system('genometk nucleotide --silent %s %s' %
                  (genome_file, output_dir))
        os.system('genometk gene --silent %s %s %s' %
                  (genome_file, gff_file, output_dir))
        os.system('genometk ssu --silent %s %s %s %s' % (genome_file,
                                                         ConfigMetadata.GTDB_SSU_GG_DB,
                                                         ConfigMetadata.GTDB_SSU_GG_TAXONOMY,
                                                         os.path.join(output_dir, ConfigMetadata.GTDB_SSU_GG_OUTPUT_DIR)))
        os.system('genometk ssu --silent %s %s %s %s' % (genome_file,
                                                         ConfigMetadata.GTDB_SSU_SILVA_DB,
                                                         ConfigMetadata.GTDB_SSU_SILVA_TAXONOMY,
                                                         os.path.join(output_dir, ConfigMetadata.GTDB_SSU_SILVA_OUTPUT_DIR)))

    def storeMetadata(self, genome_dir, genome_id=None, cur=None):
        """Parse metadata files for genome and store in database.
        Parameters
        ----------
        genome_dir : str
            Directory containing metadata files to parse.
        """

        #cur = self.conn.cursor()
        #genome_id = os.path.basename(os.path.normpath(genome_dir))
        cur.execute(
            "INSERT INTO metadata_nucleotide (id) VALUES ({0})".format(genome_id))
        cur.execute(
            "INSERT INTO metadata_genes (id) VALUES ({0})".format(genome_id))
        cur.execute(
            "INSERT INTO metadata_taxonomy (id) VALUES ({0})".format(genome_id))

        # nucleotide metadata
        metadata_nt_path = os.path.join(
            genome_dir, ConfigMetadata.GTDB_NT_FILE)
        genome_list_nt = [tuple(line.rstrip().split('\t'))
                          for line in open(metadata_nt_path)]

        query_nt = "UPDATE metadata_nucleotide SET %s = %s WHERE id = {0}".format(
            genome_id)
        print query_nt
        cur.executemany(query_nt, [(AsIs(c), v) for (c, v) in genome_list_nt])

        try:
            # protein metadata
            metadata_gene_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_GENE_FILE)
            genome_list_gene = [tuple(line.rstrip().split('\t'))
                                for line in open(metadata_gene_path)]
            query_gene = "UPDATE metadata_genes SET %s = %s WHERE id = {0}".format(
                genome_id)
            cur.executemany(query_gene, [(AsIs(c), v)
                                         for (c, v) in genome_list_gene])

            # Greengenes SSU metadata
            metadata_ssu_gg_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_SSU_GG_OUTPUT_DIR, ConfigMetadata.GTDB_SSU_FILE)
            genome_list_taxonomy, ssu_count = self._parse_taxonomy_file(
                metadata_ssu_gg_path, ConfigMetadata.GTDB_SSU_GG_PREFIX)
            query_taxonomy = "UPDATE metadata_taxonomy SET %s = %s WHERE id = {0}".format(
                genome_id)
            cur.executemany(
                query_taxonomy, [(AsIs(c), v) for (c, v) in genome_list_taxonomy])

            # SILVA SSU metadata
            metadata_ssu_silva_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_SSU_SILVA_OUTPUT_DIR, ConfigMetadata.GTDB_SSU_FILE)
            genome_list_taxonomy, ssu_count = self._parse_taxonomy_file(
                metadata_ssu_silva_path, ConfigMetadata.GTDB_SSU_SILVA_PREFIX)
            cur.executemany(
                query_taxonomy, [(AsIs(c), v) for (c, v) in genome_list_taxonomy])
            query_gene_ssu = "UPDATE metadata_genes SET ssu_count = %s WHERE id = {0}".format(
                genome_id)
            cur.execute(query_gene_ssu, (ssu_count,))
        except psycopg2.Error as e:
            raise GenomeDatabaseError(e.pgerror)

        # self.conn.commit()

    def _parse_taxonomy_file(self, metadata_taxonomy_file, prefix):
        """Parse metadata file with taxonomic information for 16S rRNA genes.

        Parameters
        ----------
        metadata_taxonomy_file : str
          Full path to file containing 16S rRNA metadata.
        Prefix : str
          Prefix to append to metadata fields.

        Returns
        -------
        list of tuples
          Field, value pairs for all metadata items.
        int
            Number of 16S sequences.
        """

        if not os.path.exists(metadata_taxonomy_file):
            return None, 0

        metadata = []
        with open(metadata_taxonomy_file) as f:
            header_line = f.readline().rstrip()
            headers = [
                prefix + '_' + x.replace('ssu_', '') for x in header_line.split('\t')]

            # Report hit to longest 16S rRNA gene. It is possible that
            # the HMMs identified a putative 16S rRNA gene, but that
            # there was no valid BLAST hit.
            longest_query_len = 0
            longest_ssu_hit_info = None
            ssu_count = 0
            for line in f:
                ssu_count += 1
                line_split = line.strip().split('\t')
                query_len = int(line_split[2])
                if query_len > longest_query_len:
                    longest_query_len = query_len
                    longest_ssu_hit_info = line_split

            if longest_ssu_hit_info:
                metadata = [(headers[i], value)
                            for i, value in enumerate(line_split)]

        return metadata, ssu_count
