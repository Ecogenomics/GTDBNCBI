###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

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
            raise e

    def exportMetadata(self, path):
        '''
        Function: exportMetadata
        Export metadata for all genomes to a csv file

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
            raise e

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
            cur.execute(query, (genome_id, meta_value))
            self.conn.commit()
            cur.close()
            self.conn.ClosePostgresConnection()
        except GenomeDatabaseError as e:
            raise e
        except psycopg2.Error as e:
            raise GenomeDatabaseError(e.pgerror)

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
                    if all_col_dict.get(key) == value['table']:
                        query_comment = "COMMENT ON COLUMN {0}.{1} IS '{2}'".format(
                            value['table'], key, value['desc'])
                        cur.execute(query_comment)
                    else:
                        logging.warning("Column {0} is already presents in the {1} table .".format(
                            key, all_col_dict.get(key)))
                else:
                    query_add_col = "ALTER TABLE {0} ADD COLUMN {1} {2}".format(
                        value['table'], key, value['type'])
                    cur.execute(query_add_col)
                    query_add_comment = "COMMENT ON COLUMN {0}.{1} IS '{2}'".format(
                        value['table'], key, value['desc'])
                    cur.execute(query_add_comment)

# ---------- PSQL are not refresh automatically so we need to drop the existing view and recreate it with a new Definition.
            cur.execute("SELECT refreshView()")
            self.conn.commit()
            cur.close()
            self.conn.ClosePostgresConnection
        except GenomeDatabaseError as e:
            raise e

    def addMetadata(self, cur, db_genome_id, genome_file, gff_file, checkm_results, output_dir):
        """Calculate and add metadata to DB.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        db_genome_id : str
            Unique database identifer of genome.
        genome_file : str
            Name of FASTA file containing nucleotide sequences.
        gff_file : str
            Name of generic feature file describing genes.
        checkm_results : dict
            CheckM metadata.
        output_dir : str
            Output directory.
        """

        # create rows for genome in metadata tables
        cur.execute(
            "INSERT INTO metadata_nucleotide (id) VALUES ({0})".format(db_genome_id))
        cur.execute(
            "INSERT INTO metadata_genes (id) VALUES ({0})".format(db_genome_id))
        cur.execute(
            "INSERT INTO metadata_taxonomy (id) VALUES ({0})".format(db_genome_id))

        self._calculateMetadata(genome_file, gff_file, output_dir)
        self._storeMetadata(cur, db_genome_id, output_dir)
        self._storeCheckM(cur, db_genome_id, checkm_results)

    def _calculateMetadata(self, genome_file, gff_file, output_dir):
        """Calculate metadata for new genome.

        Parameters
        ----------
        genome_file : str
            Name of FASTA file containing nucleotide sequences.
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

    def _storeMetadata(self, cur, db_genome_id, genome_dir):
        """Parse metadata files for genome and store in database.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        db_genome_id : str
            Unique database identifer of genome.
        genome_dir : str
            Directory containing metadata files to parse.
        """

        # nucleotide metadata
        metadata_nt_path = os.path.join(
            genome_dir, ConfigMetadata.GTDB_NT_FILE)
        genome_list_nt = [tuple(line.rstrip().split('\t'))
                          for line in open(metadata_nt_path)]
        query_nt = "UPDATE metadata_nucleotide SET %s = %s WHERE id = {0}".format(
            db_genome_id)
        cur.executemany(query_nt, [(AsIs(c), v) for (c, v) in genome_list_nt])

        try:
            # protein metadata
            metadata_gene_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_GENE_FILE)
            genome_list_gene = [tuple(line.rstrip().split('\t'))
                                for line in open(metadata_gene_path)]
            query_gene = "UPDATE metadata_genes SET %s = %s WHERE id = {0}".format(
                db_genome_id)
            cur.executemany(query_gene, [(AsIs(c), v)
                                         for (c, v) in genome_list_gene])

            # Greengenes SSU metadata
            metadata_ssu_gg_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_SSU_GG_OUTPUT_DIR, ConfigMetadata.GTDB_SSU_FILE)
            genome_list_taxonomy, ssu_count = self._parse_taxonomy_file(
                metadata_ssu_gg_path, ConfigMetadata.GTDB_SSU_GG_PREFIX)
            query_taxonomy = "UPDATE metadata_taxonomy SET %s = %s WHERE id = {0}".format(
                db_genome_id)
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
                db_genome_id)
            cur.execute(query_gene_ssu, (ssu_count,))
        except psycopg2.Error as e:
            raise GenomeDatabaseError(e.pgerror)

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

    def _storeCheckM(self, cur, db_genome_id, checkm_results):
        """Store CheckM results for genome.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        db_genome_id : str
            Unique database identifer of genome.
        checkm_results : dict
            CheckM metadata.
        """

        checkm_data = [('checkm_completeness', checkm_results['completeness']),
                       ('checkm_contamination',
                        checkm_results['contamination']),
                       ('checkm_marker_count', checkm_results['marker_count']),
                       ('checkm_marker_lineage', checkm_results['lineage']),
                       ('checkm_genome_count', checkm_results['genome_count']),
                       ('checkm_marker_set_count',
                        checkm_results['set_count']),
                       ('checkm_strain_heterogeneity', checkm_results['heterogeneity'])]

        query = "UPDATE metadata_genes SET %s = %s WHERE id = {0}".format(
            db_genome_id)
        cur.executemany(query, [(AsIs(c), v) for (c, v) in checkm_data])
