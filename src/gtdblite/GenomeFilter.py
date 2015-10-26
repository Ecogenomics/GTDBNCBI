import os
import sys
import logging

import psycopg2 as pg

from gtdblite import Tools


class GenomeFilter(object):
    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def FilterTreeData(self, GenomeDatabase, marker_ids, genome_ids,
                       comp_threshold, cont_threshold,
                       taxa_filter,
                       guaranteed_genome_list_ids, guaranteed_genome_ids,
                       individual,
                       directory, prefix):
        """Filter genomes based on provided criteria.

        Parameters
        ----------

        Returns
        -------
        set
            Database identifiers of retained genomes.
        """

        self.logger.info('Filtering genomes.')

        if not os.path.isdir(directory):
            GenomeDatabase.ReportError("Directory doesn't exist: " + directory)
            return None

        cur = GenomeDatabase.conn.cursor()

        # For all of the markers, get the expected marker size.
        if GenomeDatabase.debugMode:
            self.logger.info("Getting expected marker sizes...")

        cur.execute("SELECT markers.id, markers.name, description, id_in_database, size, external_id_prefix " +
                    "FROM markers, marker_databases " +
                    "WHERE markers.id in %s "
                    "AND markers.marker_database_id = marker_databases.id "
                    "ORDER by external_id_prefix ASC, id_in_database ASC", (tuple(marker_ids),))

        chosen_markers = dict()
        chosen_markers_order = []

        for marker_id, marker_name, marker_description, id_in_database, size, external_id_prefix in cur:
            chosen_markers[marker_id] = {'external_id_prefix': external_id_prefix, 'name': marker_name, 'description': marker_description, 'id_in_database': id_in_database, 'size': size}
            chosen_markers_order.append(marker_id)

        # output the marker info and multi_hit info
        marker_info_fh = open(os.path.join(directory, prefix + "_markers.info"), 'wb')
        multi_hits_fh = open(os.path.join(directory, prefix + "_multi_hits.info"), 'wb')
        multi_hits_header = ["Genome_ID"]
        for marker_id in chosen_markers_order:
            external_id = chosen_markers[marker_id]['external_id_prefix'] + "_" + chosen_markers[marker_id]['id_in_database']
            out_str = "\t".join([
                external_id,
                chosen_markers[marker_id]['name'],
                chosen_markers[marker_id]['description'],
                str(chosen_markers[marker_id]['size'])
            ]) + "\n"
            marker_info_fh.write(out_str)
            multi_hits_header.append(external_id)
        marker_info_fh.close()
        multi_hits_fh.write("\t".join(multi_hits_header) + "\n")

        # find genomes that are in the guaranteed list
        self.logger.info('Identified genomes to be excluded from filtering.')
        guaranteed_genomes = set()
        if guaranteed_genome_ids:
            guaranteed_genome_ids = [x.strip() for x in guaranteed_genome_ids.split(",")]
            guaranteed_genomes.update(GenomeDatabase.ExternalGenomeIdsToGenomeIds(guaranteed_genome_ids))

        if guaranteed_genome_list_ids:
            guaranteed_genome_list_ids = [x.strip() for x in guaranteed_genome_list_ids.split(",")]
            guaranteed_genomes.update(GenomeDatabase.GetGenomeIdListFromGenomeListIds(guaranteed_genome_list_ids))
        self.logger.info('Identified %d genomes to be excluded from filtering.' % len(guaranteed_genomes))

        # find genomes that fail completeness and contamination thresholds
        self.logger.info('Filtering genomes with completeness <%.1f%% or contamination >%.1f%%.' % (comp_threshold, cont_threshold))

        cur.execute("SELECT id " +
                     "FROM metadata_genes " +
                     "WHERE id IN %s " +
                     "AND (checkm_completeness < %s " +
                          "OR checkm_contamination > %s)",
                    (tuple(genome_ids), comp_threshold, cont_threshold))

        filtered_genomes = set([x[0] for x in cur])
        filtered_genomes -= guaranteed_genomes
        filtered_genomes

        self.logger.info('Filtered %d genomes based on completeness and contamination.' % len(filtered_genomes))

        if GenomeDatabase.debugMode:
            self.logger.info("Outputting aligned markers...\n")

        # select all genomes of interest (including those to be filtered)
        cur.execute("SELECT * " +
                    "FROM metadata_view "
                    "WHERE id IN %s", (tuple(genome_ids),))
        col_headers = [desc[0] for desc in cur.description]
        metadata = cur.fetchall()

        # create ARB import filter
        arb_import_filter = os.path.join(directory, prefix + "_arb_filter.ift")
        self._arbImportFilter(col_headers[2:], arb_import_filter)

        # run through each of the genomes and make the magic happen
        self.logger.info('Writing concatenated alignment and genome metadata.')
        arb_metadata_fh = open(os.path.join(directory, prefix + "_arb_metadata.txt"), 'wb')
        fasta_concat_fh = open(os.path.join(directory, prefix + "_concatenated.faa"), 'wb')
        individual_marker_fasta = dict()
        for tup in metadata:
            db_genome_id = tup[0]
            external_genome_id = tup[1]

            # filter out genomes
            if db_genome_id in filtered_genomes:
                if GenomeDatabase.debugMode:
                    checkm_completeness = tup[col_headers.index('checkm_completeness')]
                    checkm_contamination = tup[col_headers.index('checkm_contamination')]

                    self.logger.info("Genome id %s excluded. Failed CheckM quality filtering. Genomes completeness and contamination: (%f, %f)\n" %
                                     (external_genome_id, checkm_completeness, checkm_contamination))
                continue

            # write out genome info
            genome_info = dict()
            genome_info['markers'] = dict()
            genome_info['multiple_hits'] = dict()

            cur.execute("SELECT aligned_markers.marker_id, sequence, multiple_hits " +
                        "FROM aligned_markers " +
                        "WHERE genome_id = %s " +
                        "AND sequence is NOT NULL " +
                        "AND marker_id in %s ", (db_genome_id, tuple(marker_ids,)))

            if (cur.rowcount == 0):
                self.logger.warning("Genome %s has no markers for this marker set in the database and will be missing from the output files.\n" % external_id)
                continue

            for marker_id, sequence, multiple_hits in cur:
                genome_info['markers'][marker_id] = sequence
                genome_info['multiple_hits'][marker_id] = multiple_hits

            aligned_seq = ''
            multi_hits_details = []
            for marker_id in chosen_markers_order:
                multiple_hits = genome_info['multiple_hits'][marker_id]
                multi_hits_details.append('Multiple' if multiple_hits else 'Single')

                if (marker_id in genome_info['markers']) and not multiple_hits:
                    sequence = genome_info['markers'][marker_id]
                    fasta_outstr = ">%s\n%s\n" % (external_genome_id, sequence)

                    try:
                        individual_marker_fasta[marker_id].append(fasta_outstr)
                    except KeyError:
                        individual_marker_fasta[marker_id] = [fasta_outstr]
                else:
                    sequence = chosen_markers[marker_id]['size'] * '-'
                aligned_seq += sequence

            fasta_outstr = ">%s\n%s\n" % (external_genome_id, aligned_seq)
            fasta_concat_fh.write(fasta_outstr)

            multi_hits_outstr = external_genome_id + '\t' + '\t'.join(multi_hits_details) + '\n'
            multi_hits_fh.write(multi_hits_outstr)

            # write out ARB record
            multiple_hit_count = sum([1 if x == "Multiple" else 0 for x in multi_hits_details])
            self._arbRecord(arb_metadata_fh, external_genome_id,
                                col_headers[2:], tup[2:],
                                multiple_hit_count, len(chosen_markers_order),
                                aligned_seq)

        arb_metadata_fh.close()
        fasta_concat_fh.close()
        multi_hits_fh.close()

        if individual:
            self.logger.info('Writing individual alignments.')
            for marker_id in chosen_markers.keys():
                fasta_individual_fh = open(os.path.join(directory, prefix + "_" + chosen_markers[marker_id]['id_in_database'] + ".faa"), 'wb')
                fasta_individual_fh.write(''.join(individual_marker_fasta[marker_id]))
                fasta_individual_fh.close()

        return (set(genome_ids) - filtered_genomes)

    def _arbImportFilter(self, metadata_fields, output_file):
        """Create ARB import filter.

        Parameters
        ----------
        metadata_fields : list
            Names of fields to import.
        output_file : str
            Name of output file.
        """
        fout = open(output_file, 'w')
        fout.write('AUTODETECT\t"BEGIN"\n\n')
        fout.write('BEGIN\t"BEGIN*"\n\n')

        for field in ['db_name'] + metadata_fields:
            fout.write('MATCH\t"%s\\=*"\n' % field)
            fout.write('\tSRT "*\\=="\n')
            fout.write('\tWRITE "%s"\n\n' % field)

        fout.write('SEQUENCEAFTER\t"warning"\nSEQUENCESRT\t"*\\=="\nSEQUENCEEND\t"END"\n\nEND\t"END"\n')

        fout.close()

    def _arbRecord(self, fout, external_genome_id,
                        metadata_fields, metadata_values,
                        multiple_hit_count, num_marker_genes,
                        aligned_seq):
        """Write out ARB record for genome."""

        fout.write("BEGIN\n")
        fout.write("db_name=%s\n" % external_genome_id)
        for col_header, value in zip(metadata_fields, metadata_values):
            fout.write("%s=%s\n" % (col_header, str(value)))
        fout.write("multiple_homologs=%i/%i" % (multiple_hit_count, num_marker_genes))
        fout.write("aligned_seq=%s" % (aligned_seq))
        fout.write("END\n\n")
