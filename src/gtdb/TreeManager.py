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
import sys
import logging
import random
from collections import Counter, defaultdict

from numpy import (mean as np_mean,
                   std as np_std)

from biolib.taxonomy import Taxonomy

from gtdb.GenomeManager import GenomeManager
from gtdb.GenomeListManager import GenomeListManager
from gtdb.Exceptions import GenomeDatabaseError


class TreeManager(object):
    """Manages genomes, concatenated alignment, and metadata for tree inference and visualization."""

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

    def _taxa_filter(self, taxa_filter, genome_ids, guaranteed_ids, retain_guaranteed):
        """Filter genomes to specified taxa."""

        self.logger.info(
            'Filtering genomes outside taxonomic groups of interest (%s).' % taxa_filter)
        taxa_to_retain = [x.strip() for x in taxa_filter.split(',')]
        genome_ids_from_taxa = self._genomesFromTaxa(
            genome_ids, taxa_to_retain)

        retained_guaranteed_ids = set(guaranteed_ids) - genome_ids_from_taxa
        if retain_guaranteed:
            if len(retained_guaranteed_ids):
                self.logger.warning('Retaining %d guaranteed genomes from taxa not specified by the taxa filter.' % len(
                    retained_guaranteed_ids))
                self.logger.warning(
                    "You can use the '--guaranteed_taxa_filter' flag to filter these genomes.")

            genomes_to_retain = genome_ids.intersection(
                genome_ids_from_taxa).union(guaranteed_ids)
        else:
            genomes_to_retain = genome_ids.intersection(genome_ids_from_taxa)
            self.logger.info("Filtered %d 'guaranteed' genomes based on taxonomic affiliations." % len(
                retained_guaranteed_ids))

        self.logger.info('Filtered %d genomes based on taxonomic affiliations.' % (
            len(genome_ids) - len(genomes_to_retain)))

        return genomes_to_retain

    def filterGenomes(self, marker_ids,
                      genome_ids,
                      quality_threshold,
                      quality_weight,
                      comp_threshold,
                      cont_threshold,
                      min_perc_aa,
                      min_rep_perc_aa,
                      taxa_filter,
                      guaranteed_taxa_filter,
                      genomes_to_exclude,
                      guaranteed_ids,
                      rep_ids,
                      directory,
                      prefix):
        """Filter genomes based on provided criteria.

        Parameters
        ----------

        Returns
        -------
        set
            Database identifiers of retained genomes.
        """

        if not os.path.exists(directory):
            os.makedirs(directory)

        # get mapping from db genome IDs to external IDs
        genome_mngr = GenomeManager(self.cur, self.currentUser)
        external_ids = genome_mngr.genomeIdsToExternalGenomeIds(genome_ids)
        filter_genome_file = os.path.join(
            directory, prefix + '_filtered_genomes.tsv')
        fout_filtered = open(filter_genome_file, 'w')

        self.logger.info('Filtering initial set of %d genomes.' %
                         len(genome_ids))

        extra_guaranteed_ids = [
            x for x in guaranteed_ids if x not in genome_ids]
        if len(extra_guaranteed_ids) > 0:
            self.logger.warning('Identified {0} guaranteed genomes absent from specified input genomes (Those genomes will not appear in the final tree).'.format(
                len(extra_guaranteed_ids)))
            guaranteed_ids = [x for x in guaranteed_ids if x in genome_ids]
        self.logger.info(
            'Identified %d genomes to be excluded from filtering.' % len(guaranteed_ids))

        # for all markers, get the expected marker size
        self.cur.execute("SELECT markers.id, markers.name, description, id_in_database, size, external_id_prefix " +
                         "FROM markers, marker_databases " +
                         "WHERE markers.id in %s "
                         "AND markers.marker_database_id = marker_databases.id "
                         "ORDER by external_id_prefix ASC, id_in_database ASC", (tuple(marker_ids),))

        chosen_markers = dict()
        chosen_markers_order = []

        total_alignment_len = 0
        for marker_id, marker_name, marker_description, id_in_database, size, external_id_prefix in self.cur:
            chosen_markers[marker_id] = {'external_id_prefix': external_id_prefix, 'name': marker_name,
                                         'description': marker_description, 'id_in_database': id_in_database, 'size': size}
            chosen_markers_order.append(marker_id)
            total_alignment_len += size

        # filter genomes based on taxonomy
        genomes_to_retain = genome_ids
        if taxa_filter:
            new_genomes_to_retain = self._taxa_filter(taxa_filter,
                                                      genomes_to_retain,
                                                      guaranteed_ids,
                                                      retain_guaranteed=True)
            for genome_id in genomes_to_retain - new_genomes_to_retain:
                rep_str = 'Representative' if genome_id in rep_ids else ''
                fout_filtered.write('%s\t%s\t%s\n' % (
                    external_ids[genome_id], 'Filtered on taxonomic affiliation.', rep_str))

            genomes_to_retain = new_genomes_to_retain

        if guaranteed_taxa_filter:
            new_genomes_to_retain = self._taxa_filter(guaranteed_taxa_filter,
                                                      genomes_to_retain,
                                                      guaranteed_ids,
                                                      retain_guaranteed=False)
            for genome_id in genomes_to_retain - new_genomes_to_retain:
                rep_str = 'Representative' if genome_id in rep_ids else ''
                fout_filtered.write('%s\t%s\t%s\n' % (
                    external_ids[genome_id], 'Filtered on guaranteed taxonomic affiliation.', rep_str))

            genomes_to_retain = new_genomes_to_retain

        # find genomes based on completeness, contamination, or genome quality
        self.logger.info('Filtering genomes with completeness <%.1f%%, contamination >%.1f%%, or quality <%.1f%% (weight = %.1f).' % (
            comp_threshold,
            cont_threshold,
            quality_threshold,
            quality_weight))
        filtered_genomes = self._filterOnGenomeQuality(genomes_to_retain,
                                                       quality_threshold,
                                                       quality_weight,
                                                       comp_threshold,
                                                       cont_threshold)

        # sanity check representatives are not of poor quality
        final_filtered_genomes = set()
        for genome_id, quality in filtered_genomes.iteritems():
            if genome_id not in guaranteed_ids:
                if genome_id in rep_ids:
                    self.logger.warning('Retaining representative genome %s despite poor estimated quality (comp=%.1f%%, cont=%.1f%%).' % (
                        external_ids[genome_id], quality[0], quality[1]))
                else:
                    final_filtered_genomes.add(genome_id)
                    fout_filtered.write(
                        '%s\t%s\t%.2f\t%.2f\n' % (external_ids[genome_id],
                                                  'Filtered on quality (completeness, contamination).',
                                                  quality[0],
                                                  quality[1]))

        self.logger.info('Filtered %d genomes based on completeness, contamination, and quality.' % len(
            final_filtered_genomes))

        genomes_to_retain -= final_filtered_genomes

        # filter genomes explicitly specified for exclusion
        if genomes_to_exclude:
            for genome_id in genomes_to_exclude:
                if genome_id in external_ids:
                    fout_filtered.write('%s\t%s\n' % (
                        external_ids[genome_id], 'Explicitly marked for exclusion.'))

            conflicting_genomes = guaranteed_ids.intersection(
                genomes_to_exclude)
            if conflicting_genomes:
                raise GenomeDatabaseError('Genomes marked for both retention and exclusion, e.g.: %s'
                                          % conflicting_genomes.pop())

            new_genomes_to_retain = genomes_to_retain.difference(
                genomes_to_exclude)
            self.logger.info('Filtered %d genomes explicitly indicated for exclusion.' % (
                len(genomes_to_retain) - len(new_genomes_to_retain)))
            genomes_to_retain = new_genomes_to_retain

        # filter genomes with insufficient number of amino acids in MSA
        self.logger.info(
            'Filtering genomes with insufficient amino acids in the MSA.')
        filter_on_aa = set()
        for genome_id in genomes_to_retain:
            aligned_marker_query = ("SELECT sequence, multiple_hits,hit_number,unique_genes " +
                                    "FROM aligned_markers " +
                                    "WHERE genome_id = %s " +
                                    "AND sequence is NOT NULL " +
                                    "AND marker_id IN %s")

            self.cur.execute(aligned_marker_query,
                             (genome_id, tuple(marker_ids)))

            total_aa = 0
            for sequence, multiple_hits, hit_number, unique_genes in self.cur:
                if not multiple_hits:
                    total_aa += len(sequence) - sequence.count('-')
                elif unique_genes == 1:
                    total_aa += len(sequence) - sequence.count('-')

            # should retain guaranteed genomes unless they have zero amino
            # acids in MSA
            if genome_id in guaranteed_ids:
                if total_aa != 0:
                    continue
                else:
                    self.logger.warning(
                        'Filtered guaranteed genome %s with zero amino acids in MSA.' % external_ids[genome_id])

            perc_alignment = total_aa * 100.0 / total_alignment_len
            if perc_alignment < min_perc_aa:
                rep_str = ''
                if genome_id in rep_ids:
                    if perc_alignment < min_rep_perc_aa:
                        rep_str = 'Representative'
                        self.logger.warning('Filtered representative genome %s due to lack of aligned amino acids (%.1f%%).' % (
                            external_ids[genome_id], perc_alignment))
                    else:
                        self.logger.warning('Retaining representative genome %s despite small numbers of aligned amino acids (%.1f%%).' % (
                            external_ids[genome_id], perc_alignment))
                        continue

                filter_on_aa.add(genome_id)
                fout_filtered.write(
                    '%s\t%s\t%d\t%.1f\t%s\n' % (external_ids[genome_id],
                                                'Insufficient number of amino acids in MSA (total AA, % alignment length)',
                                                total_aa,
                                                perc_alignment,
                                                rep_str))

        fout_filtered.close()

        self.logger.info(
            'Filtered %d genomes with insufficient amino acids in the MSA.' % len(filter_on_aa))

        genomes_to_retain.difference_update(filter_on_aa)
        self.logger.info('Producing tree data for %d genomes.' %
                         len(genomes_to_retain))

        good_genomes_file = os.path.join(
            directory, prefix + '_good_genomes.tsv')
        good_genomes = open(good_genomes_file, 'w')
        for item in genomes_to_retain:
            good_genomes.write("{0}\n".format(item))
        good_genomes.close()

        return (genomes_to_retain, chosen_markers_order, chosen_markers)

    def _mimagQualityInfo(self, metadata, col_headers):
        """Add MIMAG quality information to metadata."""

        col_headers.append('mimag_high_quality')
        col_headers.append('mimag_medium_quality')
        col_headers.append('mimag_low_quality')

        comp_index = col_headers.index('checkm_completeness')
        cont_index = col_headers.index('checkm_contamination')

        gtdb_domain_index = col_headers.index('gtdb_domain')
        lsu_5s_length_index = col_headers.index('lsu_5s_length')
        lsu_23s_length_index = col_headers.index('lsu_silva_length')
        ssu_length_index = col_headers.index('ssu_silva_length')
        trna_aa_count_index = col_headers.index('trna_aa_count')

        for i, gm in enumerate(metadata):
            comp = float(gm[comp_index])
            cont = float(gm[cont_index])

            lsu_5s_length = 0
            if gm[lsu_5s_length_index]:
                lsu_5s_length = int(gm[lsu_5s_length_index])

            lsu_23s_length = 0
            if gm[lsu_23s_length_index]:
                lsu_23s_length = int(gm[lsu_23s_length_index])

            ssu_length = 0
            if gm[ssu_length_index]:
                ssu_length = int(gm[ssu_length_index])

            trna_aa_count = 0
            if gm[trna_aa_count_index]:
                trna_aa_count = int(gm[trna_aa_count_index])

            gtdb_domain = gm[gtdb_domain_index]
            ssu_length_threshold = 1200
            if gtdb_domain == 'd__Archaea':
                ssu_length_threshold = 900

            hq = False
            mq = False
            lq = False
            if comp > 90 and cont < 5:
                if (ssu_length >= ssu_length_threshold
                        and lsu_23s_length >= 1900
                        and lsu_5s_length >= 80
                        and trna_aa_count_index >= 18):
                    hq = True
                else:
                    mq = True
            elif comp >= 50 and cont < 10:
                mq = True
            elif cont < 10:
                lq = True

            gm += (hq, mq, lq)
            metadata[i] = gm

        return metadata

    def writeFiles(self,
                   marker_ids,
                   genomes_to_retain,
                   cols_per_gene,
                   min_consensus,
                   max_consensus,
                   rnd_seed,
                   min_perc_taxa,
                   min_perc_aa,
                   chosen_markers_order,
                   chosen_markers,
                   alignment,
                   individual,
                   directory,
                   prefix,
                   no_trim,
                   classic_header):
        '''
        Write summary files and ARB files

        :param marker_ids:
        :param genomes_to_retain:
        :param chosen_markers_order:
        :param chosen_markers:
        :param alignment:
        :param individual:
        :param directory:
        :param prefix:
        '''

        if not os.path.exists(directory):
            os.makedirs(directory)

        # output the marker info and multiple hit info
        multi_hits_fh = open(
            os.path.join(directory, prefix + "_multi_hits.tsv"), 'wb')
        multi_hits_header = ["Genome_ID"]
        for marker_id in chosen_markers_order:
            external_id = chosen_markers[marker_id][
                'external_id_prefix'] + "_" + chosen_markers[marker_id]['id_in_database']
            multi_hits_header.append(external_id)
        multi_hits_fh.write("\t".join(multi_hits_header) + "\n")

        # select genomes to retain
        self.cur.execute("SELECT * " +
                         "FROM metadata_view "
                         "WHERE id IN %s", (tuple(genomes_to_retain),))
        col_headers = [desc[0] for desc in self.cur.description]
        metadata = self.cur.fetchall()

        # add MIMAG quality information
        # metadata = self._mimagQualityInfo(metadata, col_headers)

        # identify columns of interest
        genome_accn_idx = col_headers.index('accession')
        # to export MSA with original header ( not canonical)
        if classic_header:
            genome_name_index = col_headers.index('accession')
        else:
            genome_name_index = col_headers.index('formatted_accession')
        genome_id_index = col_headers.index('id')
        col_headers.remove('id')
        # col_headers.remove('formatted_accession')

        # run through each of the genomes and concatenate markers
        self.logger.info(
            'Concatenated marker genes for %d genomes.' % len(genomes_to_retain))

        individual_marker_fasta = dict()
        single_copy = defaultdict(int)
        multi_copy_identical = defaultdict(int)
        ubiquitous = defaultdict(int)
        multi_hits_details = defaultdict(list)
        msa = {}

        for genome_metadata in metadata:
            db_genome_id = genome_metadata[genome_id_index]
            external_genome_id = genome_metadata[genome_name_index]

            # get aligned markers
            aligned_marker_query = ("SELECT am.marker_id, sequence, multiple_hits, evalue, unique_genes " +
                                    "FROM aligned_markers am " +
                                    "LEFT JOIN markers m on m.id=am.marker_id " +
                                    "WHERE genome_id = %s " +
                                    "AND sequence is NOT NULL " +
                                    "AND marker_id in %s " +
                                    "ORDER BY m.id_in_database")

            self.cur.execute(aligned_marker_query,
                             (db_genome_id, tuple(marker_ids)))

            if (self.cur.rowcount == 0):
                self.logger.warning(
                    "Genome %s has no markers for this marker set and will be missing from the output files." % external_genome_id)
                continue

            genome_info = dict()
            genome_info['markers'] = dict()
            genome_info['multiple_hits'] = dict()
            genome_info['unique_genes'] = dict()
            for marker_id, sequence, multiple_hits, evalue, unique_genes in self.cur:
                if evalue:  # markers without an e-value are missing
                    genome_info['markers'][marker_id] = sequence
                genome_info['multiple_hits'][marker_id] = multiple_hits
                if unique_genes is not None:
                    genome_info['unique_genes'][marker_id] = int(unique_genes)
                else:
                    genome_info['unique_genes'][marker_id] = 0

            aligned_seq = ''
            for marker_id in chosen_markers_order:
                multiple_hits = genome_info['multiple_hits'][marker_id]
                unique_genes = genome_info['unique_genes'][marker_id]

                if (marker_id in genome_info['markers']):
                    ubiquitous[marker_id] += 1
                    if not multiple_hits:
                        single_copy[marker_id] += 1
                        multi_hits_details[db_genome_id].append('Single')
                    elif multiple_hits and unique_genes == 1:
                        multi_copy_identical[marker_id] += 1
                        multi_hits_details[db_genome_id].append(
                            'Multi_Identical')
                    else:
                        multi_hits_details[db_genome_id].append('Multiple')
                else:
                    multi_hits_details[db_genome_id].append('Missing')

                if (marker_id in genome_info['markers']) and (not multiple_hits or (multiple_hits and unique_genes == 1)):
                    sequence = genome_info['markers'][marker_id]
                    fasta_outstr = ">%s\n%s\n" % (external_genome_id, sequence)

                    try:
                        individual_marker_fasta[marker_id].append(fasta_outstr)
                    except KeyError:
                        individual_marker_fasta[marker_id] = [fasta_outstr]
                else:
                    sequence = chosen_markers[marker_id]['size'] * '-'
                    fasta_outstr = ">%s\n%s\n" % (external_genome_id, sequence)

                    try:
                        individual_marker_fasta[marker_id].append(fasta_outstr)
                    except KeyError:
                        individual_marker_fasta[marker_id] = [fasta_outstr]
                aligned_seq += sequence

            msa[external_genome_id] = aligned_seq
            multi_hits_outstr = '%s\t%s\n' % (
                external_genome_id, '\t'.join(multi_hits_details[db_genome_id]))
            multi_hits_fh.write(multi_hits_outstr)

        multi_hits_fh.close()

        # filter columns without sufficient representation across taxa
        if no_trim:
            self.logger.info('Trimming step is skipped.')
            trimmed_seqs = msa

        else:
            self.logger.info(
                'Trimming columns with insufficient taxa or poor consensus.')
            mask, pruned_seqs = self._trim_seqs(cols_per_gene,
                                                min_perc_aa / 100.0,
                                                min_consensus / 100.0,
                                                max_consensus / 100.0,
                                                min_perc_taxa / 100.0,
                                                rnd_seed,
                                                msa,
                                                chosen_markers_order,
                                                chosen_markers)

            trimmed_seqs = dict(pruned_seqs)

            # write out mask for MSA
            msa_mask_out = open(os.path.join(
                directory, prefix + "_mask.txt"), 'w')
            msa_mask_out.write(''.join(['1' if m else '0' for m in mask]))
            msa_mask_out.close()

        # write out MSA
        fasta_concat_filename = os.path.join(
            directory, prefix + "_concatenated.faa")
        fasta_concat_fh = open(fasta_concat_filename, 'wb')
        for genome_id, aligned_seq in trimmed_seqs.iteritems():
            fasta_outstr = ">%s\n%s\n" % (genome_id, aligned_seq)
            fasta_concat_fh.write(fasta_outstr)
        fasta_concat_fh.close()

        # create dictionary indicating NCBI species and clustered genomes
        self.logger.info(
            'Determining NCBI species assignments for genomes in GTDB clusters.')
        self.cur.execute("SELECT accession, ncbi_taxonomy, gtdb_clustered_genomes " +
                         "FROM metadata_view")

        ncbi_sp = {}
        gtdb_clustered_gids = {}
        for row in self.cur.fetchall():
            gid = row[0]
            ncbi_sp[gid] = row[1].split(';')[-1].strip()

            if row[2] is not None:
                gtdb_clustered_gids[gid] = row[2].split(';')
            else:
                gtdb_clustered_gids[gid] = []

        # write out ARB metadata and ARB import filter
        self.logger.info(
            'Writing ARB metadata for %d genomes.' % len(genomes_to_retain))
        arb_metadata_file = os.path.join(
            directory, prefix + "_arb_metadata.txt")
        arb_metadata_fh = open(arb_metadata_file, 'wb')

        col_headers.append('gtdb_cluster_ncbi_species')

        arb_import_filter = os.path.join(directory, prefix + "_arb_filter.ift")
        self._arbImportFilter(col_headers, arb_import_filter)

        for genome_metadata in metadata:
            # take special care of the genome identifier and name as these
            # are handle as a special case in the ARB metadata file
            genome_metadata = list(genome_metadata)

            db_genome_id = genome_metadata[genome_id_index]
            genome_accn = genome_metadata[genome_accn_idx]
            external_genome_id = genome_metadata[genome_name_index]
            del genome_metadata[genome_id_index]
            # ==================================================================
            # del genome_metadata[max(genome_id_index, genome_name_index)]
            # print(genome_metadata[0:10])
            # del genome_metadata[min(genome_id_index, genome_name_index)]
            # print(genome_metadata[0:10])
            # ==================================================================

            # add NCBI species information field to ARB record
            ncbi_sp_count = defaultdict(int)
            for gid in gtdb_clustered_gids[genome_accn]:
                ncbi_sp_count[ncbi_sp[gid]] += 1

            ncbi_sp_str = []
            for sp, count in sorted(ncbi_sp_count.items(), key=lambda kv: kv[1], reverse=True):
                ncbi_sp_str.append('%s: %d' % (sp, count))
            genome_metadata.append('; '.join(ncbi_sp_str))

            # write out ARB record
            multiple_hit_count = sum(
                [1 if x == "Multiple" else 0 for x in multi_hits_details[db_genome_id]])
            msa_gene_count = sum(
                [1 if x == "Single" else 0 for x in multi_hits_details[db_genome_id]])

            aligned_seq = trimmed_seqs[external_genome_id]
            if not alignment:
                aligned_seq = ''
            self._arbRecord(arb_metadata_fh,
                            external_genome_id,
                            col_headers,
                            genome_metadata,
                            multiple_hit_count,
                            msa_gene_count,
                            len(chosen_markers_order),
                            aligned_seq)

        arb_metadata_fh.close()

        # write out marker gene summary info
        marker_info_fh = open(
            os.path.join(directory, prefix + "_markers_info.tsv"), 'wb')
        marker_info_fh.write(
            'Marker Id\tName\tDescription\tLength (bp)\tSingle copy (%)\tUbiquity (%)\n')

        for marker_id in chosen_markers_order:
            external_id = chosen_markers[marker_id][
                'external_id_prefix'] + "_" + chosen_markers[marker_id]['id_in_database']

            sc = single_copy.get(marker_id, 0) * 100.0 / len(genomes_to_retain)
            u = ubiquitous.get(marker_id, 0) * 100.0 / len(genomes_to_retain)

            out_str = "\t".join([
                external_id,
                chosen_markers[marker_id]['name'],
                chosen_markers[marker_id]['description'],
                str(chosen_markers[marker_id]['size']),
                '%.2f' % sc,
                '%.2f' % u
            ]) + "\n"
            marker_info_fh.write(out_str)
        marker_info_fh.close()

        # write out individual marker gene alignments
        if individual:
            self.logger.info('Writing individual alignments.')
            for marker_id in chosen_markers.keys():
                fasta_individual_fh = open(os.path.join(
                    directory, prefix + "_" + chosen_markers[marker_id]['id_in_database'] + ".faa"), 'wb')
                fasta_individual_fh.write(
                    ''.join(individual_marker_fasta[marker_id]))
                fasta_individual_fh.close()

        return fasta_concat_filename

    def _trim_seqs(self, cols_per_gene, min_perc_aa,
                   min_consensus,
                   max_consensus,
                   min_per_taxa,
                   rnd_seed,
                   msa,
                   chosen_markers_order,
                   chosen_markers):
        """Randomly select a subset of columns from the MSA of each marker."""

        markers_sizeinfo = []

        for chosen_marker in chosen_markers_order:
            markers_sizeinfo.append((chosen_markers.get(
                chosen_marker).get('id_in_database'),
                chosen_markers.get(
                chosen_marker).get('name'),
                chosen_markers.get(chosen_marker).get('size')))

        max_gaps = 1.0 - min_per_taxa

        # randomly select columns meeting filtering criteria
        self.logger.info(
            'Randomly sampling %d columns passing filtering criteria from each marker gene.' % cols_per_gene)

        mask, output_seqs = self.subsample_msa(
            msa, markers_sizeinfo, cols_per_gene, max_gaps, min_consensus, max_consensus, rnd_seed)
        return mask, output_seqs

    def subsample_msa(self, seqs, markers, cols_per_gene, max_gaps, min_identical_aa, max_identical_aa, rnd_seed):
        """Sample columns from each marker in multiple sequence alignment."""

        alignment_length = len(seqs.values()[0])
        sampled_cols = []
        start = 0
        lack_sufficient_cols = 0
        lack_cols_marker_ids = []
        avg_perc_cols = []

        count_wrong_pa = 0
        count_wrong_cons = 0

        random.seed(rnd_seed)

        for marker_id, marker_name, marker_len in markers:
            end = start + marker_len

            valid_cols, count_wrong_pa, count_wrong_cons = self.identify_valid_columns(marker_name, count_wrong_pa,
                                                                                       count_wrong_cons,
                                                                                       start,
                                                                                       end,
                                                                                       seqs,
                                                                                       max_gaps,
                                                                                       min_identical_aa,
                                                                                       max_identical_aa)

            assert(len(valid_cols) <= marker_len)  # sanity check

            self.logger.info('%s: S:%d, E:%d, LEN:%d, COLS:%d, PERC:%.1f' % (
                marker_name,
                start,
                end,
                marker_len,
                len(valid_cols),
                len(valid_cols) * 100.0 / marker_len))

            avg_perc_cols.append(len(valid_cols) * 100.0 / marker_len)

            if len(valid_cols) < cols_per_gene:
                self.logger.warning(
                    'Marker has <%d columns after filtering.' % cols_per_gene)
                lack_sufficient_cols += 1
                lack_cols_marker_ids.append(marker_id)

            offset_valid_cols = [i + start for i in valid_cols]
            sampled_cols.extend(random.sample(
                offset_valid_cols, min(cols_per_gene, len(offset_valid_cols))))

            start = end
        mask = [1 if i in sampled_cols else 0 for i in range(alignment_length)]

        self.logger.info('Identified %d of %d marker genes with <%d columns for sampling:' % (
            lack_sufficient_cols,
            len(markers),
            cols_per_gene))
        self.logger.info('%s' % ', '.join(lack_cols_marker_ids))
        self.logger.info('Marker genes had %.1f+/-%.1f%% of columns available for selection on average.' % (
            np_mean(avg_perc_cols),
            np_std(avg_perc_cols)))

        # trim columns
        output_seqs = {}
        for seq_id, seq in seqs.iteritems():
            masked_seq = ''.join([seq[i]
                                  for i in xrange(0, len(mask)) if mask[i]])
            output_seqs[seq_id] = masked_seq

        self.logger.info('Trimmed alignment from %d to %d AA (%d by minimum taxa percent, %d by consensus, maximum of %d columns per genes).' % (len(seqs[seqs.keys()[0]]),
                                                                                                                                                 len(output_seqs[output_seqs.keys()[0]]), count_wrong_pa, count_wrong_cons, cols_per_gene))
        self.logger.info('Final MSA contains %d columns.' % len(sampled_cols))

        return mask, output_seqs

    def identify_valid_columns(self, name, count_wrong_pa, count_wrong_cons, start, end, seqs, max_gaps, min_identical_aa, max_identical_aa):
        """Identify columns meeting gap and amino acid ubiquity criteria."""

        GAP_CHARS = set(['-', '.', '_', '*'])
        STANDARD_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')

        gap_count = defaultdict(int)
        amino_acids = [list() for _ in xrange(end - start)]
        num_genomes = 0
        for seq_id, seq in seqs.iteritems():
            num_genomes += 1
            gene_seq = seq[start:end].upper()
            for i, ch in enumerate(gene_seq):
                if ch in GAP_CHARS:
                    gap_count[i] += 1
                else:
                    amino_acids[i].append(ch)

        valid_cols = set()
        for i in xrange(0, end - start):
            if float(gap_count.get(i, 0)) / num_genomes <= max_gaps:
                c = Counter(amino_acids[i])

                if not c.most_common(1):
                    continue

                letter, count = c.most_common(1)[0]
                if letter not in STANDARD_AMINO_ACIDS:
                    self.logger.warning(
                        'Most common amino acid was not in standard alphabet: %s' % letter)

                aa_ratio = float(count) / (num_genomes - gap_count.get(i, 0))
                if min_identical_aa <= aa_ratio < max_identical_aa:
                    valid_cols.add(i)
                else:
                    count_wrong_cons += 1
            else:
                count_wrong_pa += 1

        return valid_cols, count_wrong_pa, count_wrong_cons

    def _filterOnGenomeQuality(self, genome_ids, quality_threshold, quality_weight, comp_threshold, cont_threshold):
        """Filter genomes on completeness and contamination thresholds.

        Parameters
        ----------
        genome_ids : list
            Database identifier for genomes of interest.
        quality_threshold : float
            Minimum required quality threshold.
        quality_weight : float
            Weighting factor for assessing genome quality.
        comp_threshold : float
            Minimum required completeness.
        cont_threshold : float
            Maximum permitted contamination.

        Returns
        -------
        set
            Database identifier of genomes failing quality filtering.
        """

        self.cur.execute("SELECT id, checkm_completeness, checkm_contamination " +
                         "FROM metadata_genes " +
                         "WHERE id IN %s " +
                         "AND (checkm_completeness < %s " +
                         "OR checkm_contamination > %s " +
                         "OR (checkm_completeness - %s*checkm_contamination) < %s)",
                         (tuple(genome_ids), comp_threshold, cont_threshold, quality_weight, quality_threshold))

        return {x[0]: [x[1], x[2]] for x in self.cur}

    def _genomesFromTaxa(self, genome_ids, taxa_to_retain):
        """Filter genomes to those within specified taxonomic groups.

        Parameters
        ----------
        genome_ids : list
            Database identifier for genomes of interest.
        taxa_to_retain : list
            Taxonomic groups of interest.

        Returns
        -------
        set
            Database identifier of genomes from specified taxonomic groups.
        """

        taxa_to_retain_at_rank = [[]
                                  for _ in xrange(len(Taxonomy.rank_prefixes))]
        for taxon in taxa_to_retain:
            taxon_prefix = taxon[0:3]
            if taxon_prefix not in Taxonomy.rank_prefixes:
                self.logger.error('Invalid taxon prefix: %s' % taxon)
                sys.exit(0)

            rank_index = Taxonomy.rank_index[taxon_prefix]
            taxa_to_retain_at_rank[rank_index].append(taxon)

        query_str = []
        query_tuple = [tuple(genome_ids)]
        for i, taxa in enumerate(taxa_to_retain_at_rank):
            if len(taxa):
                query_str.append('gtdb_%s IN %%s' % Taxonomy.rank_labels[i])
                query_tuple.append(tuple(taxa))
        query_str = ' OR '.join(query_str)

        self.cur.execute("SELECT id " +
                         "FROM metadata_taxonomy " +
                         "WHERE id IN %s "
                         "AND " + query_str,
                         tuple(query_tuple))

        genome_ids_from_taxa = set([x[0] for x in self.cur])

        return genome_ids_from_taxa

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

        fout.write('MATCH\t"%s\\=*"\n' % 'db_name')
        fout.write('\tSRT "*\\=="\n')
        fout.write('\tWRITE "%s"\n\n' % 'name')

        # place organism name near top for convenience
        fout.write('MATCH\t"%s\\=*"\n' % 'organism_name')
        fout.write('\tSRT "*\\=="\n')
        fout.write('\tWRITE "%s"\n\n' % 'organism_name')

        ncbi_taxonomy_index = metadata_fields.index('ncbi_taxonomy')
        metadata_fields.insert(ncbi_taxonomy_index + 1, 'ncbi_genus')
        metadata_fields.insert(ncbi_taxonomy_index + 2, 'ncbi_species')
        fields = metadata_fields + ['msa_gene_count',
                                    'msa_num_marker_genes',
                                    'msa_aa_count',
                                    'msa_length',
                                    'multiple_homologs']
        for field in fields:
            if field != 'organism_name':
                fout.write('MATCH\t"%s\\=*"\n' % field)
                fout.write('\tSRT "*\\=="\n')
                fout.write('\tWRITE "%s"\n\n' % field)

        fout.write('SEQUENCEAFTER\t"multiple_homologs*"\n')
        fout.write('SEQUENCESRT\t"*\\=="\n')
        fout.write('SEQUENCEEND\t"END"\n\n')
        fout.write('END\t"END"\n')

        fout.close()

    def canonical_gid(self, gid):
        """Get canonical form of NCBI genome accession.

        Example:
            G005435135 -> G005435135
            GCF_005435135.1 -> G005435135
            GCF_005435135.1_ASM543513v1_genomic -> G005435135
            RS_GCF_005435135.1 -> G005435135
            GB_GCA_005435135.1 -> G005435135
        """

        if gid.startswith('U'):
            return gid

        gid = gid.replace('RS_', '').replace('GB_', '')
        gid = gid.replace('GCA_', 'G').replace('GCF_', 'G')
        if '.' in gid:
            gid = gid[0:gid.find('.')]

        return gid

    def _arbRecord(self, fout,
                   external_genome_id,
                   metadata_fields,
                   metadata_values,
                   multiple_hit_count,
                   msa_gene_count,
                   num_marker_genes,
                   aligned_seq,
                   classic_header=False):
        """Write out ARB record for genome."""

        # customize output relative to raw database table
        if (classic_header and (external_genome_id.startswith('GB') or external_genome_id.startswith('RS'))) or ( not classic_header and external_genome_id.startswith('G')):
            metadata_values = list(metadata_values)
            metadata_fields = list(metadata_fields)
            ncbi_taxonomy_index = metadata_fields.index('ncbi_taxonomy')
            ncbi_genus = metadata_values[ncbi_taxonomy_index].split(';')[
                5].strip()
            ncbi_species = metadata_values[ncbi_taxonomy_index].split(';')[
                6].strip()
            if 'ncbi_genus' not in metadata_fields:
                metadata_fields.insert(ncbi_taxonomy_index + 1, 'ncbi_genus')
                metadata_values.insert(ncbi_taxonomy_index + 1, ncbi_genus)
            else:
                metadata_values.insert(
                    metadata_fields.index('ncbi_genus'), ncbi_genus)

            if 'ncbi_species' not in metadata_fields:
                metadata_fields.insert(ncbi_taxonomy_index + 2, 'ncbi_species')
                metadata_values.insert(ncbi_taxonomy_index + 2, ncbi_species)
            else:
                metadata_values.insert(metadata_fields.index(
                    'ncbi_species'), ncbi_species)
            organism_name_index = metadata_fields.index('organism_name')
            ncbi_organism_name_index = metadata_fields.index(
                'ncbi_organism_name')
            metadata_values[organism_name_index] = metadata_values[ncbi_organism_name_index]

            gtdb_clustered_genomes_index = metadata_fields.index(
                'gtdb_clustered_genomes')
            if metadata_values[gtdb_clustered_genomes_index] != None:
                metadata_values[gtdb_clustered_genomes_index] = ';'.join([self.canonical_gid(
                    x) for x in metadata_values[gtdb_clustered_genomes_index].split(';')])

        fout.write("BEGIN\n")
        fout.write("db_name=%s\n" % external_genome_id)
        for col_header, value in zip(metadata_fields, metadata_values):
            if isinstance(value, float):
                value = '%.4g' % value

            # replace equal signs as these are incompatible with the ARB parser
            value = str(value)
            if value:
                value = value.replace('=', '/')

            fout.write("%s=%s\n" % (col_header, value))

        fout.write("msa_gene_count=%d\n" % msa_gene_count)
        fout.write("msa_num_marker_genes=%d\n" % num_marker_genes)
        fout.write("msa_aa_count=%d\n" %
                   (len(aligned_seq) - aligned_seq.count('-')))
        fout.write("msa_length=%d\n" % len(aligned_seq))
        fout.write("multiple_homologs=%d\n" % multiple_hit_count)
        fout.write("aligned_seq=%s\n" % (aligned_seq))
        fout.write("END\n\n")
