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
import sys

import psycopg2
from psycopg2.extensions import AsIs

import ConfigMetadata
import Config
from Exceptions import GenomeDatabaseError

from biolib.taxonomy import Taxonomy
from biolib.seq_io import read_fasta


class MetadataManager(object):

    def __init__(self, cur, currentUser):
        """Initialize.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        currentUser : User
            Current user of database.
        """

        self.logger = logging.getLogger()

        self.cur = cur
        self.currentUser = currentUser

    def viewMetadata(self):
        '''
         Function: viewMetadata
        Lists all metadata field available in the database

        Returns:
        Print lists of Metadata on Screen
        '''
        try:
            self.cur.execute("SELECT * FROM view_list_meta_columns")
            print "\t".join(("Table", "Field", "Datatype", "Description"))
            for tup in self.cur.fetchall():
                info_list = list(tup)
                if info_list[2] == "double precision":
                    info_list[2] = "float"
                elif str(info_list[2]).startswith("timestamp"):
                    info_list[2] = "timestamp"
                print "\t".join(info_list)
        except GenomeDatabaseError as e:
            raise e

    def exportMetadata(self, path):
        '''
        Function: exportMetadata
        Export metadata for all genomes to a csv file

        :param path: Path to the output file
        '''

        try:
            query_tmp = "CREATE TEMPORARY TABLE temp_tb as SELECT * FROM metadata_view;"
            query_tmp += "ALTER TABLE temp_tb DROP COLUMN id;ALTER TABLE temp_tb DROP COLUMN study_id;"
            self.cur.execute(query_tmp)
            query = "SELECT * FROM temp_tb"
            outputquery = 'copy ({0}) to stdout with csv header'.format(query)
            with open(path, 'w') as f:
                self.cur.copy_expert(outputquery, f)
            print "Export Successful"
        except GenomeDatabaseError as e:
            raise e

    def ExportGenomePaths(self, path):
        '''
        Function: exportMetadata
        Export the full path for all genomes to a csv file

        :param path: Path to the output file
        '''

        try:
            query_tmp = "SELECT id_at_source,external_id_prefix,fasta_file_location FROM genomes,genome_sources WHERE genomes.genome_source_id=genome_sources.id;"
            self.cur.execute(query_tmp)
            with open(path, "w") as f:
                for (id, prefix, file_location) in self.cur:
                    dir_prefix = None
                    if prefix == 'U':
                        dir_prefix = Config.GTDB_GENOME_USR_DIR
                    elif prefix == 'RS':
                        dir_prefix = Config.GTDB_GENOME_RSQ_DIR
                    elif prefix == 'GB':
                        dir_prefix = Config.GTDB_GENOME_GBK_DIR
                    else:
                        raise GenomeDatabaseError(
                            "Unrecognized database prefix: %s" % prefix)
                    f.write("{0}\t{1}\n".format(prefix + "_" + id, os.path.dirname(os.path.join(dir_prefix, file_location))))
            print "Export Successful"
        except GenomeDatabaseError as e:
            raise e

    def exportTaxonomy(self, taxonomy_src, output_file):
        """Write taxonomy to file.

        Parameters
        ----------
        taxonomy_src : str
          Indicates desired taxonomy ('GTDB' or 'NCBI').
        output_file : str
          Output file.
        """

        try:
            if taxonomy_src == 'NCBI':
                self.cur.execute("SELECT genome, ncbi_taxonomy FROM metadata_view")
            else:
                self.cur.execute("SELECT genome, gtdb_taxonomy FROM metadata_view")

            fout = open(output_file, 'w')
            for genome_id, taxonomy in self.cur.fetchall():
                if taxonomy:
                    fout.write('%s\t%s\n' % (genome_id, taxonomy))
                else:
                    fout.write('%s\t%s\n' % (genome_id, ';'.join(Taxonomy.rank_prefixes)))
            fout.close()

            print 'Taxonomy information written to: %s' % output_file
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
            self.cur.execute(query, (genome_id, meta_value))
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
            self.cur.execute(query)
            all_col_dict = dict(self.cur.fetchall())
            for key, value in data_dict.iteritems():
                if key in all_col_dict:
                    if all_col_dict.get(key) == value['table']:
                        query_comment = "COMMENT ON COLUMN {0}.{1} IS '{2}'".format(
                            value['table'], key, value['desc'])
                        self.cur.execute(query_comment)
                    else:
                        logging.warning("Column {0} is already presents in the {1} table .".format(
                            key, all_col_dict.get(key)))
                else:
                    query_add_col = "ALTER TABLE {0} ADD COLUMN {1} {2}".format(
                        value['table'], key, value['type'])
                    self.cur.execute(query_add_col)
                    query_add_comment = "COMMENT ON COLUMN {0}.{1} IS '{2}'".format(
                        value['table'], key, value['desc'])
                    self.cur.execute(query_add_comment)

# ---------- PSQL are not refresh automatically so we need to drop the existing view and recreate it with a new Definition.
            self.cur.execute("SELECT refreshView()")

        except GenomeDatabaseError as e:
            raise e

    def addMetadata(self, db_genome_id, genome_file, gff_file, checkm_results, output_dir):
        """Calculate and add metadata to DB.

        Parameters
        ----------
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
        self.cur.execute(
            "INSERT INTO metadata_nucleotide (id) VALUES ({0})".format(db_genome_id))
        self.cur.execute(
            "INSERT INTO metadata_genes (id) VALUES ({0})".format(db_genome_id))
        self.cur.execute(
            "INSERT INTO metadata_taxonomy (id) VALUES ({0})".format(db_genome_id))
        self.cur.execute(
            "INSERT INTO metadata_rna (id) VALUES ({0})".format(db_genome_id))
        self.cur.execute(
            "INSERT INTO metadata_sequence (id) VALUES ({0})".format(db_genome_id))

        self._calculateMetadata(genome_file, gff_file, output_dir)
        self._storeMetadata(db_genome_id, output_dir)
        self._storeCheckM(db_genome_id, checkm_results)
        return True

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

        os.system('genometk rna --silent --db %s --taxonomy_file %s %s ssu %s' % (
            ConfigMetadata.GTDB_SSU_GG_DB,
            ConfigMetadata.GTDB_SSU_GG_TAXONOMY,
            genome_file,
            os.path.join(output_dir, ConfigMetadata.GTDB_SSU_GG_OUTPUT_DIR)))

        os.system('genometk rna --silent --db %s --taxonomy_file %s %s ssu %s' % (
            ConfigMetadata.GTDB_SSU_SILVA_DB,
            ConfigMetadata.GTDB_SSU_SILVA_TAXONOMY,
            genome_file,
            os.path.join(output_dir, ConfigMetadata.GTDB_SSU_SILVA_OUTPUT_DIR)))

        os.system('genometk rna --silent --db %s --taxonomy_file %s %s lsu_23S %s' % (
            ConfigMetadata.GTDB_LSU_SILVA_DB,
            ConfigMetadata.GTDB_LSU_SILVA_TAXONOMY,
            genome_file,
            os.path.join(output_dir, ConfigMetadata.GTDB_LSU_SILVA_OUTPUT_DIR)))
        return True

    def _storeMetadata(self, db_genome_id, genome_dir):
        """Parse metadata files for genome and store in database.

        Parameters
        ----------
        db_genome_id : str
            Unique database identifier of genome.
        genome_dir : str
            Directory containing metadata files to parse.
        """
        try:
            # nucleotide metadata
            metadata_nt_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_NT_FILE)
            genome_list_nt = [tuple(line.rstrip().split('\t'))
                              for line in open(metadata_nt_path)]
            query_nt = "UPDATE metadata_nucleotide SET %s = %s WHERE id = {0}".format(
                db_genome_id)
            for c, v in genome_list_nt:
                try:
                    v = float(v)
                    self.cur.execute(query_nt, [AsIs(c), v])
                except:
                    self.cur.execute(query_nt, [AsIs(c), v])

            # protein metadata
            metadata_gene_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_GENE_FILE)
            genome_list_gene = [tuple(line.rstrip().split('\t'))
                                for line in open(metadata_gene_path)]
            query_gene = "UPDATE metadata_genes SET %s = %s WHERE id = {0}".format(
                db_genome_id)
            for c, v in genome_list_gene:
                try:
                    v = float(v)
                    self.cur.execute(query_gene, [AsIs(c), v])
                except:
                    self.cur.execute(query_gene, [AsIs(c), v])

            # Greengenes SSU metadata
            query_taxonomy = "UPDATE metadata_rna SET %s = %s WHERE id = {0}".format(
                db_genome_id)
            query_sequence = "UPDATE metadata_sequence SET %s = %s WHERE id = {0}".format(
                db_genome_id)
            metadata_ssu_gg_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_SSU_GG_OUTPUT_DIR, ConfigMetadata.GTDB_SSU_FILE)
            genome_list_taxonomy, _ssu_count, ssu_query_id = self._parse_taxonomy_file(
                metadata_ssu_gg_path, ConfigMetadata.GTDB_SSU_GG_PREFIX)
            if genome_list_taxonomy:
                for c, v in genome_list_taxonomy:
                    try:
                        if "blast_subject_id" not in c:
                            v = float(v)
                        self.cur.execute(query_taxonomy, [AsIs(c), v])
                    except:
                        self.cur.execute(query_taxonomy, [AsIs(c), v])

            # SILVA SSU metadata saved in metadata_ssu table [HACK: eventually information will only be stored in this table]
            query_taxonomy = "UPDATE metadata_rna SET %s = %s WHERE id = {0}".format(
                db_genome_id)
            query_sequence = "UPDATE metadata_sequence SET %s = %s WHERE id = {0}".format(
                db_genome_id)
            metadata_ssu_silva_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_SSU_SILVA_OUTPUT_DIR, ConfigMetadata.GTDB_SSU_FILE)
            metadata_ssu_fna_silva_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_SSU_SILVA_OUTPUT_DIR, ConfigMetadata.GTDB_SSU_FNA_FILE)
            metadata_ssu_silva_summary_file = os.path.join(genome_dir, ConfigMetadata.GTDB_SSU_SILVA_OUTPUT_DIR, ConfigMetadata.GTDB_SSU_SILVA_SUMMARY_FILE)

            genome_list_taxonomy, ssu_count, ssu_query_id = self._parse_taxonomy_file(
                metadata_ssu_silva_path, ConfigMetadata.GTDB_SSU_SILVA_PREFIX, metadata_ssu_silva_summary_file)
            if genome_list_taxonomy:
                for c, v in genome_list_taxonomy:
                    try:
                        if "blast_subject_id" not in c:
                            v = float(v)
                        self.cur.execute(query_taxonomy, [AsIs(c), v])
                    except:
                        self.cur.execute(query_taxonomy, [AsIs(c), v])
                if ssu_query_id is not None:
                    genome_list_sequence = self._parse_sequence_file(metadata_ssu_fna_silva_path, ConfigMetadata.GTDB_SSU_SILVA_PREFIX, ssu_query_id)
                    for c, v in genome_list_sequence:
                        self.cur.execute(query_sequence, [AsIs(c), v])

            # SILVA LSU metadata saved in metadata_ssu table [HACK: eventually information will only be stored in this table]
            query_taxonomy = "UPDATE metadata_rna SET %s = %s WHERE id = {0}".format(
                db_genome_id)
            metadata_lsu_silva_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_LSU_SILVA_OUTPUT_DIR, ConfigMetadata.GTDB_LSU_FILE)
            metadata_lsu_fna_silva_path = os.path.join(
                genome_dir, ConfigMetadata.GTDB_LSU_SILVA_OUTPUT_DIR, ConfigMetadata.GTDB_LSU_FNA_FILE)
            metadata_lsu_silva_summary_file = os.path.join(genome_dir, ConfigMetadata.GTDB_LSU_SILVA_OUTPUT_DIR, ConfigMetadata.GTDB_LSU_SILVA_SUMMARY_FILE)
            genome_list_taxonomy, lsu_count, lsu_query_id = self._parse_taxonomy_file(
                metadata_lsu_silva_path, ConfigMetadata.GTDB_LSU_SILVA_PREFIX, metadata_lsu_silva_summary_file)
            if genome_list_taxonomy:
                for c, v in genome_list_taxonomy:
                    try:
                        if "blast_subject_id" not in c:
                            v = float(v)
                        self.cur.execute(query_taxonomy, [AsIs(c), v])
                    except:
                        self.cur.execute(query_taxonomy, [AsIs(c), v])
                if lsu_query_id is not None:
                    genome_list_sequence = self._parse_sequence_file(metadata_lsu_fna_silva_path, ConfigMetadata.GTDB_LSU_SILVA_PREFIX, lsu_query_id)
                    for c, v in genome_list_sequence:
                        self.cur.execute(query_sequence, [AsIs(c), v])

            query_gene_ssu = "UPDATE metadata_genes SET ssu_count = %s WHERE id = {0}".format(db_genome_id)
            self.cur.execute(query_gene_ssu, (ssu_count,))

            query_gene_lsu = "UPDATE metadata_genes SET lsu_count = %s WHERE id = {0}".format(db_genome_id)
            self.cur.execute(query_gene_lsu, (lsu_count,))

            return True
        except psycopg2.Error as e:
            print "error"
            raise GenomeDatabaseError(e.pgerror)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

    def _parse_taxonomy_file(self, metadata_taxonomy_file, prefix, summary_file=None):
        """Parse metadata file with taxonomic information for rRNA genes.

        Parameters
        ----------
        metadata_taxonomy_file : str
          Full path to file containing 16S rRNA metadata.
        Prefix : str
          Prefix to append to metadata fields.
        Fna file : str
          Full path to file containing 16S rRNA sequences.

        Returns
        -------
        list of tuples
          Field, value pairs for all metadata items.
        int
            Number of 16S sequences.
        """

        if not os.path.exists(metadata_taxonomy_file):
            return None, 0, None

        metadata = []
        ssu_query_id = None

        num_lines = sum(1 for line in open(metadata_taxonomy_file))
        if num_lines == 1:  # Empty file
            return metadata, 0, None

        with open(metadata_taxonomy_file) as f:
            header_line = f.readline().rstrip()
            headers = [prefix + '_' + x for x in header_line.split('\t')]

            split_headers = header_line.rstrip().split("\t")
            for pos in range(0, len(split_headers)):
                header = split_headers[pos]
                if header == 'query_id':
                    query_id_pos = pos
                    break

            # Report hit to longest 16S rRNA gene. It is possible that
            # the HMMs identified a putative 16S rRNA gene, but that
            # there was no valid BLAST hit.
            longest_query_len = 0
            longest_ssu_hit_info = None
            identified_ssu_genes = 0
            for line in f:
                line_split = line.strip().split('\t')
                query_len = int(line_split[2])
                if query_len > longest_query_len:
                    longest_query_len = query_len
                    longest_ssu_hit_info = line_split
                    ssu_query_id = line_split[query_id_pos]

            if longest_ssu_hit_info:
                metadata = [(headers[i], value)
                            for i, value in enumerate(longest_ssu_hit_info)]

                if summary_file is not None and os.path.exists(summary_file):
                    with open(summary_file) as fsum:
                        header_line = fsum.readline()  # consume header line
                        header_list = [x.strip() for x in header_line.split('\t')]
                        idx_seq = header_list.index("Sequence length")
                        for line in fsum:
                            identified_ssu_genes += 1
                            sum_list = [x.strip() for x in line.split('\t')]
                            if sum_list[0] == ssu_query_id:
                                metadata.append(('{0}_contig_len'.format(prefix), sum_list[idx_seq]))

            return metadata, identified_ssu_genes, ssu_query_id

    def _parse_sequence_file(self, fna_file, prefix, ssu_query_id):
        metadata = []
        all_genes_dict = read_fasta(fna_file, False)
        sequence = all_genes_dict[ssu_query_id]
        metadata.append(('{0}_sequence'.format(prefix), sequence))
        return metadata

    def _storeCheckM(self, db_genome_id, checkm_results):
        """Store CheckM results for genome.

        Parameters
        ----------
        db_genome_id : str
            Unique database identifer of genome.
        checkm_results : dict
            CheckM metadata.
        """
        try:
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
            self.cur.executemany(query, [(AsIs(c), v) for (c, v) in checkm_data])
            return True
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise
