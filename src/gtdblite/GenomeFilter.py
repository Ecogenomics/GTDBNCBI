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

        gg_fh = open(os.path.join(directory, prefix + "_concatenated.arbtxt"), 'wb')
        fasta_concat_fh = open(os.path.join(directory, prefix + "_concatenated.faa"), 'wb')
        marker_info_fh = open(os.path.join(directory, prefix + "_concatenated_markers.info"), 'wb')
        multi_hits_fh = open(os.path.join(directory, prefix + "_multi_hits.info"), 'wb')

        # output the marker info and multi_hit info
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

        cur.execute("SELECT external_id_prefix || '_' || id_at_source as external_id, genomes.id, checkm_completeness, checkm_contamination " +
                     "FROM genomes " +
                     "LEFT OUTER JOIN genome_sources ON genomes.genome_source_id = genome_sources.id " +
                     "WHERE genomes.id in %s " +
                     "AND (checkm_completeness < %s " +
                          "OR checkm_contamination > %s)",
                    (tuple(genome_ids), comp_threshold, cont_threshold))

        filtered_genomes = set([values[1] for values in cur])
        filtered_genomes -= guaranteed_genomes

        self.logger.info('Filtered %d genomes based on completeness and contamination.' % len(filtered_genomes))

        if GenomeDatabase.debugMode:
            self.logger.info("Outputting aligned markers...\n")

        # Select all the genomes - could be lots of genome id, create a temp table
        temp_table_name = Tools.generateTempTableName()

        cur.execute("CREATE TEMP TABLE %s (id integer)" % (temp_table_name,))
        query = "INSERT INTO {0} (id) VALUES (%s)".format(temp_table_name)
        cur.executemany(query, [(x,) for x in genome_ids])

        query = ("SELECT external_id_prefix || '_' || id_at_source as external_id, genomes.id, genomes.name, checkm_completeness, checkm_contamination, owned_by_root, users.username " +
                 "FROM genomes " +
                 "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
                 "LEFT OUTER JOIN genome_sources ON genomes.genome_source_id = genome_sources.id " +
                 "WHERE genomes.id in (SELECT id from {0})".format(temp_table_name))

        cur.execute(query)

        genome_details = cur.fetchall()

        # run through each of the genomes and make the magic happen
        self.logger.info('Writing concatenated alignment and genome metadata.')
        individual_marker_fasta = dict()
        genome_info = dict()
        for (external_id, genome_id, name, checkm_completeness, checkm_contamination, owned_by_root, owner) in genome_details:
            if genome_id in filtered_genomes:
                if GenomeDatabase.debugMode:
                    self.logger.info("Genome id %s excluded. Failed CheckM quality filtering. Genomes completeness and contamination: (%f, %f)\n" %
                                     (external_id, checkm_completeness, checkm_contamination))
                continue

            # Populate genome info
            genome_info['markers'] = dict()
            genome_info['multiple_hits'] = dict()
            genome_info['name'] = name
            genome_info['external_id'] = external_id
            genome_info['owner'] = ('root' if owned_by_root else owner)

            cur.execute("SELECT aligned_markers.marker_id, sequence, multiple_hits " +
                        "FROM aligned_markers " +
                        "WHERE genome_id = %s " +
                        "AND sequence is NOT NULL " +
                        "AND marker_id in %s ", (genome_id, tuple(marker_ids,)))

            if (cur.rowcount == 0):
                sys.stderr.write("WARNING: Genome %s has no markers for this marker set in the database and will be missing from the output files.\n" % external_id)
                sys.stderr.flush()
                continue

            for marker_id, sequence, multiple_hits in cur:
                genome_info['markers'][marker_id] = sequence
                genome_info['multiple_hits'][marker_id] = multiple_hits

            aligned_seq = ''
            multi_hits_details = [genome_info['external_id']]
            for marker_id in chosen_markers_order:
                multiple_hits = genome_info['multiple_hits'][marker_id]
                multi_hits_details.append('Multiple' if multiple_hits else 'Single')

                if (marker_id in genome_info['markers']) and not multiple_hits:  # (not multiple_hits or 'include_multihits' in config_dict):
                    sequence = genome_info['markers'][marker_id]
                    fasta_outstr = ">%s\n%s\n" % (genome_info['external_id'],
                                                  sequence)
                    try:
                        individual_marker_fasta[marker_id].append(fasta_outstr)
                    except KeyError:
                        individual_marker_fasta[marker_id] = [fasta_outstr]
                else:
                    sequence = chosen_markers[marker_id]['size'] * '-'
                aligned_seq += sequence

            fasta_outstr = ">%s\n%s\n" % (genome_info['external_id'],
                                          aligned_seq)
            multi_hits_outstr = "\t".join(multi_hits_details) + "\n"

            gg_list = ["BEGIN",
                        "db_name=%s" % genome_info['external_id'],
                        "organism=%s" % genome_info['name'],
                        "prokMSA_id=%s" % genome_info['external_id'],
                        "owner=%s" % genome_info['owner'],
                        "checkm_completeness=%f" % checkm_completeness,
                        "checkm_contamination=%f" % checkm_contamination,
                        "multiple_homologs=%i/%i" % (sum([1 if x == "Multiple" else 0 for x in multi_hits_details]), len(chosen_markers_order)),  # Lazy, lazy. Fix this.
                        "warning=",
                        "aligned_seq=%s" % (aligned_seq),
                        "END"]

            gg_outstr = "\n".join(gg_list) + "\n\n"

            fasta_concat_fh.write(fasta_outstr)
            gg_fh.write(gg_outstr)
            multi_hits_fh.write(multi_hits_outstr)

        gg_fh.close()
        fasta_concat_fh.close()
        multi_hits_fh.close()

        if individual:
            self.logger.info('Writing individual alignments.')
            for marker_id in chosen_markers.keys():
                fasta_individual_fh = open(os.path.join(directory, prefix + "_" + chosen_markers[marker_id]['id_in_database'] + ".faa"), 'wb')
                fasta_individual_fh.write(''.join(individual_marker_fasta[marker_id]))
                fasta_individual_fh.close()

        return True
