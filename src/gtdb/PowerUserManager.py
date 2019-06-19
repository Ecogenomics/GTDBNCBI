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

import DefaultValues
import Config
from Exceptions import GenomeDatabaseError

from biolib.taxonomy import Taxonomy


from AlignedMarkerManager import AlignedMarkerManager
from MarkerSetManager import MarkerSetManager


class PowerUserManager(object):

    def __init__(self, cur, currentUser, release, threads=1):
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
        self.threads = threads
        self.db_release = release

    def runTreeWeightedExceptions(self, path, comp, conta, qweight, qt):
        '''
        Function: runTreeWeightedException
        Export list of NCBI genomes that do comply the filter criteria but are of interest.

        :param path: Path to the output file
        '''
        try:
            if (not self.currentUser.isRootUser()):
                raise GenomeDatabaseError(
                    "Only the root user can run this command")
                return False
            self.cur.execute("SELECT id,mt.ncbi_taxonomy FROM genomes g " +
                             "LEFT JOIN metadata_genes mg USING (id) " +
                             "LEFT JOIN metadata_taxonomy mt  USING (id) " +
                             "LEFT JOIN metadata_ncbi mn  USING (id) " +
                             "WHERE g.genome_source_id IN (2,3) and (mt.gtdb_genome_representative is not NULL or " +
                             "(mt.gtdb_genome_representative is NULL and mg.checkm_completeness > %s and mg.checkm_contamination < %s " +
                             "and mg.checkm_completeness-%s*mg.checkm_contamination > %s)) and mt.ncbi_taxonomy is not NULL",
                             (comp, conta, qweight, qt))
            print self.cur.mogrify("SELECT id,mt.ncbi_taxonomy FROM genomes g " +
                                   "LEFT JOIN metadata_genes mg USING (id) " +
                                   "LEFT JOIN metadata_taxonomy mt  USING (id) " +
                                   "LEFT JOIN metadata_ncbi mn  USING (id) " +
                                   "WHERE g.genome_source_id IN (2,3) and (mt.gtdb_genome_representative is not NULL or " +
                                   "(mt.gtdb_genome_representative is NULL and mg.checkm_completeness > %s and mg.checkm_contamination < %s " +
                                   "and mg.checkm_completeness-%s*mg.checkm_contamination > %s)) and mt.ncbi_taxonomy is not NULL",
                                   (comp, conta, qweight, qt))

            processed_results = zip(*self.cur)
            existing_id = processed_results[0]
            existing_taxonomy = processed_results[1]
            order_list = [x.split(';')[3] for x in existing_taxonomy]
            self.cur.execute("SELECT g.id,g.name,mg.checkm_completeness,mg.checkm_contamination,mt.ncbi_taxonomy,mnuc.genome_size,(mg.checkm_completeness-4*mg.checkm_contamination) as quality_threshold,mn.ncbi_organism_name " +
                             "FROM genomes g " +
                             "LEFT JOIN metadata_genes mg USING (id) " +
                             "LEFT JOIN metadata_ncbi mn  USING (id) " +
                             "LEFT JOIN metadata_nucleotide mnuc  USING (id) " +
                             "LEFT JOIN metadata_taxonomy mt  USING (id) " +
                             "WHERE g.genome_source_id IN (2,3) and " +
                             "(mg.checkm_completeness > %s and  mg.checkm_contamination < %s " +
                             "and mg.checkm_completeness-4*mg.checkm_contamination > %s) and mt.ncbi_taxonomy is not NULL and g.id not in %s",
                             (DefaultValues.EXCEPTION_FILTER_ONE_CHECKM_COMPLETENESS, DefaultValues.EXCEPTION_FILTER_ONE_CHECKM_CONTAMINATION, DefaultValues.EXCEPTION_FILTER_ONE_QUALITY_THRESHOLD, existing_id))

            dict_except_order = {}
            for (gid, name, compl, conta, ncbitax, size, qual, orga) in self.cur:
                if ncbitax.split(';')[3] != 'o__':
                    if self._checkTaxonomyUniqueness(ncbitax, order_list):
                        if self._checkTaxonomyUniqueness(ncbitax, dict_except_order):
                            dict_except_order[ncbitax.split(';')[3]] = {'quality': float(qual), 'id': gid, 'full_info': [
                                name, compl, conta, ncbitax, size, qual, orga, 'First']}
                        else:
                            if dict_except_order.get(ncbitax.split(';')[3]).get('quality') < float(qual):
                                dict_except_order[ncbitax.split(';')[3]] = {'quality': float(qual), 'id': gid, 'full_info': [
                                    name, compl, conta, ncbitax, size, qual, orga, 'First']}
                else:
                    dict_except_order[gid] = {'quality': float(qual), 'id': gid, 'full_info': [
                        name, compl, conta, ncbitax, size, qual, orga, 'unknown order']}

            fh = open(path, "w")
            fh.write(
                "Name,CheckM_Completeness,CheckM_Contamination,NCBI_Taxonomy,Genome_size,Quality_Threshold,Organism_name,Filter_passed\n")
            for _k, item in dict_except_order.iteritems():
                fh.write(",".join(str(v)
                                  for v in item.get('full_info')) + "\n")
            fh.close()

        except GenomeDatabaseError as e:
            raise e
        return True

    def _checkTaxonomyUniqueness(self, ncbitax, order_list):
        ncbi_tax_order = ncbitax.split(';')[3]
        if ncbi_tax_order in order_list:
            return False

        return True

    def _validateSSU_LSU(self):
        """Validate SSU and LSU gene metadata.

        If a gene has been identified, it should have an
        associated length. Conversely, a length should not
        be given in the gene count is zero.
        """

        query = ("SELECT accession, "
                 "ssu_count, ssu_length, "
                 "lsu_23s_count, lsu_23s_length, "
                 "lsu_5s_count, lsu_5s_length "
                 "FROM metadata_view")
        self.cur.execute(query)

        for d in self.cur:
            (gid, ssu_count, ssu_length,
                lsu_23s_count, lsu_23s_length,
                lsu_5s_count, lsu_5s_length) = d

            if ssu_count >= 1 and not ssu_length:
                print 'Missing 16S length information: %s' % gid
            if ssu_count == 0 and ssu_length:
                print 'No 16S gene identified, but length information is provided: %s' % gid

            if lsu_23s_count >= 1 and not lsu_23s_length:
                print 'Missing 23S length information: %s' % gid
            if lsu_23s_count == 0 and lsu_23s_length:
                print 'No 23S gene identified, but length information is provided: %s' % gid

            if lsu_5s_count >= 1 and not lsu_5s_length:
                print 'Missing 5S length information: %s' % gid
            if lsu_5s_count == 0 and lsu_5s_length:
                print 'No 5S gene identified, but length information is provided: %s' % gid

    def _validateTypeStrains(self):
        """Validate 'type_strain' field."""

        query = ("SELECT accession, ncbi_organism_name, type_strain, lpsn_strain, dsmz_strain, straininfo_strain "
                 "FROM metadata_view "
                 "WHERE ncbi_organism_name is not NULL")
        self.cur.execute(query)

        for gid, ncbi_organism_name, type_strain, lpsn_strain, dsmz_strain, straininfo_strain in self.cur:
            ncbi_organism_name = ncbi_organism_name.replace(
                'Candidatus ', '')     # normalize species name
            ncbi_sp = ' '.join(ncbi_organism_name.split()[
                               :2])          # get species name
            ncbi_strain_ids = set([strain_id.strip(
            ) for strain_id in ncbi_organism_name.replace(ncbi_sp, '').split('=')])
            ncbi_strain_ids_no_spaces = set()
            for ncbi_strain_id in ncbi_strain_ids:
                ncbi_strain_ids_no_spaces.add(ncbi_strain_id.replace(' ', ''))

            for authority, strains in (['lpsn', lpsn_strain], ['dsmz', dsmz_strain], ['straininfo', straininfo_strain]):
                # check cases when authority indicates a type strain
                if type_strain and authority in type_strain:
                    if not strains or ncbi_sp not in strains:
                        print 'Incorrect type strain assignment attributed to %s: %s' % (authority.upper(), gid)
                    else:
                        strain_ids = set(
                            [strain_id.strip() for strain_id in strains.replace(ncbi_sp, '').split('=')])
                        if strain_ids.intersection(ncbi_strain_ids) or strain_ids.intersection(ncbi_strain_ids_no_spaces):
                            print 'Incorrect type strain assignment attributed to %s: %s' % (authority.upper(), gid)

                # check for missing authority
                if not type_strain or authority not in type_strain:
                    if strains and ncbi_sp in strains:
                        strain_ids = set(
                            [strain_id.strip() for strain_id in strains.replace(ncbi_sp, '').split('=')])
                        if strain_ids.intersection(ncbi_strain_ids) or strain_ids.intersection(ncbi_strain_ids_no_spaces):
                            print 'Missing type strain assignment to %s: %s' % (authority.upper(), gid)

    def _validateMIMAG(self):
        """Validationg MIMAG assignments."""

        query = ("SELECT accession, mimag_high_quality, mimag_medium_quality, mimag_low_quality, "
                 "checkm_completeness, checkm_contamination, trna_aa_count, "
                 "ssu_count, ssu_length, "
                 "lsu_23s_count, lsu_23s_length, "
                 "lsu_5s_count, lsu_5s_length, "
                 "gtdb_domain "
                 "FROM metadata_view "
                 "WHERE gtdb_domain is not NULL")
        self.cur.execute(query)

        for d in self.cur:
            (gid, mimag_high_quality, mimag_medium_quality, mimag_low_quality,
                checkm_completeness, checkm_contamination, trna_aa_count,
                ssu_count, ssu_length,
                lsu_23s_count, lsu_23s_length,
                lsu_5s_count, lsu_5s_length,
                gtdb_domain) = d

            min_lsu_5s_length = 80
            min_lsu_23s_length = 1900
            if gtdb_domain == 'd__Bacteria':
                min_ssu_length = 1200
            elif gtdb_domain == 'd__Archaea':
                min_ssu_length = 900
            else:
                print 'Genome %s has an unrecognized domain assignment: %s' % (gid, gtdb_domain)

            if (checkm_completeness > 90 and checkm_contamination < 5 and
                    trna_aa_count >= 18 and
                    ssu_count >= 1 and ssu_length >= min_ssu_length and
                    lsu_23s_count >= 1 and lsu_23s_length >= min_lsu_23s_length and
                    lsu_5s_count >= 1 and lsu_5s_length >= min_lsu_5s_length):
                if not mimag_high_quality:
                    print 'Failed to mark genome %s as MIMAG high quality.' % gid
                if mimag_medium_quality:
                    print 'Incorrectly marked genome %s as MIMAG medium quality.' % gid
                if mimag_low_quality:
                    print 'Incorrectly marked genome %s as MIMAG low quality.' % gid
            elif checkm_completeness >= 50 and checkm_contamination <= 10:
                if mimag_high_quality:
                    print 'Incorrectly marked genome %s as MIMAG high quality.' % gid
                if not mimag_medium_quality:
                    print 'Failed to mark genome %s as MIMAG medium quality.' % gid
                if mimag_low_quality:
                    print 'Incorrectly marked genome %s as MIMAG low quality.' % gid
            elif checkm_contamination <= 10:
                if mimag_high_quality:
                    print 'Incorrectly marked genome %s as MIMAG high quality.' % gid
                if mimag_medium_quality:
                    print 'Incorrectly marked genome %s as MIMAG medium quality.' % gid
                if not mimag_low_quality:
                    print 'Failed to mark genome %s as MIMAG low quality.' % gid

    def runSanityCheck(self):
        try:
            if (not self.currentUser.isRootUser()):
                raise GenomeDatabaseError(
                    "Only the root user can run this command")

            # validate type strains
            self.logger.info(
                'Validating 5S, 16S, and 23S count and gene length data.')
            self._validateSSU_LSU()

            # validate type strains
            #==========We should review this later==============================
            # self.logger.info('Validating type strain.')
            # self._validateTypeStrains()
            #===================================================================

            # validate MIMAG assignments
            self.logger.info('Validating MIMAG assignments.')
            self._validateMIMAG()

            # check if the representatives are still in the database
            query = ("SELECT id FROM genomes where genome_source_id in (2,3)")
            self.cur.execute(query)
            ncbi_ids = [gid for (gid,) in self.cur]

            query = ("SELECT id,id_at_source FROM genomes")
            self.cur.execute(query)
            raw_ids = [(gid, source_id) for (gid, source_id) in self.cur]
            all_ids, all_source_ids = zip(*raw_ids)
            dict_all_ids = {k: v for (k, v) in raw_ids}

            query = (
                "SELECT distinct(gtdb_genome_representative) from metadata_taxonomy where gtdb_genome_representative is not NULL")
            self.cur.execute(query)
            representatives = [self._chompRecord(
                record) for (record,) in self.cur]

            for representative in representatives:
                if representative not in all_source_ids:
                    print "REPRESENTATIVE {0} has been removed from the database".format(representative)

            query = ("SELECT id,protein_count,ncbi_submitter from metadata_genes LEFT JOIN metadata_ncbi using (id) WHERE id in (SELECT id from genomes where genome_source_id in (2,3));")
            self.cur.execute(query)
            dict_meta_ncbi = {gid: {"count": count, "submitter": project}
                              for (gid, count, project) in self.cur}

            for ncbi_genome in ncbi_ids:
                if ncbi_genome not in dict_meta_ncbi:
                    print "{0} has no metadata in metadata_ncbi".format(dict_all_ids[ncbi_genome])
                else:
                    if dict_meta_ncbi[ncbi_genome]["count"] is None or dict_meta_ncbi[ncbi_genome]["count"] == '' or dict_meta_ncbi[ncbi_genome]["count"] == 0:
                        print "{0} protein_count value in metadata_nucleotide is {1}".format(dict_all_ids[ncbi_genome], dict_meta_ncbi[ncbi_genome]["count"])
                    if dict_meta_ncbi[ncbi_genome]["submitter"] is None or dict_meta_ncbi[ncbi_genome]["submitter"] == '':
                        print "{0} ncbi_submitter value in metadata_ncbi is {1}".format(dict_all_ids[ncbi_genome], dict_meta_ncbi[ncbi_genome]["submitter"])

            query = (
                "SELECT id,checkm_completeness,protein_count from metadata_genes")
            self.cur.execute(query)
            dict_meta_genes = {gid: {"checkm": checkm, "protein_count": count} for (
                gid, checkm, count) in self.cur}

            for genome in all_ids:
                if genome not in dict_meta_genes:
                    print "{0} has no metadata in metadata_genes".format(dict_all_ids[genome])
                else:
                    if dict_meta_genes[genome]["checkm"] is None or dict_meta_genes[genome]["checkm"] == '':
                        print "{0} checkm_completeness value in metadata_genes is {1}".format(dict_all_ids[genome], dict_meta_genes[genome]["checkm"])
                    if dict_meta_genes[genome]["protein_count"] is None or dict_meta_genes[genome]["protein_count"] == '' or dict_meta_genes[genome]["protein_count"] == 0:
                        print "{0} protein_count value in metadata_genes is {1}".format(dict_all_ids[genome], dict_meta_genes[genome]["protein_count"])

            query = ("SELECT id,gc_count from metadata_nucleotide")
            self.cur.execute(query)
            dict_meta_nuc = {gid: {"gc": gc} for (gid, gc) in self.cur}

            for genome in all_ids:
                if genome not in dict_meta_nuc:
                    print "{0} has no metadata in metadata_nucleotide".format(dict_all_ids[genome])
                else:
                    if dict_meta_nuc[genome]["gc"] is None or dict_meta_nuc[genome]["gc"] == '' or dict_meta_nuc[genome]["gc"] == 0:
                        print "{0} gc_count value in metadata_nucleotide is {1}".format(dict_all_ids[genome], dict_meta_nuc[genome]["gc"])

        except GenomeDatabaseError as e:
            raise e

        return True

    def runTaxonomyCheck(self, rank_depth):
        """Compare GTDB taxonomy to NCBI taxonomy, and report differences."""

        try:
            # Check if gtdb_domain is the same as the gtdb_taxonomy and
            # ncbi_taxonomy domain
            query = ("SELECT g.id_at_source, mt.gtdb_domain, mt.ncbi_taxonomy, gtv.gtdb_taxonomy FROM metadata_taxonomy mt " +
                     "LEFT JOIN genomes g using (id) LEFT JOIN gtdb_taxonomy_view gtv USING (id) " +
                     "WHERE genome_source_id in (2,3) ")
            self.cur.execute(query)
            list_domain = [[a, b, c, d] for (a, b, c, d) in self.cur]
            print "#Conflicting Domains:"
            for genome_id, gtdb_domain, ncbi_taxonomy, gtdb_taxonomy in list_domain:
                if gtdb_taxonomy is not None and not gtdb_taxonomy.startswith('d__;'):
                    if gtdb_domain is None or gtdb_domain == 'd__':
                        print '{0}\tgtdb_domain:{1}\tncbi_taxonomy:{2}\tgtdb_taxonomy:{3}'.format(genome_id, gtdb_domain, ncbi_taxonomy, gtdb_taxonomy)
                    elif gtdb_domain not in gtdb_taxonomy:
                        print '{0}\tgtdb_domain:{1}\tncbi_taxonomy:{2}\tgtdb_taxonomy:{3}'.format(genome_id, gtdb_domain, ncbi_taxonomy, gtdb_taxonomy)

                if ncbi_taxonomy is not None and not ncbi_taxonomy.startswith('d__;'):
                    if gtdb_domain is None or gtdb_domain == 'd__':
                        print '{0}\tgtdb_domain:{1}\tncbi_taxonomy:{2}\tgtdb_taxonomy:{3}'.format(genome_id, gtdb_domain, ncbi_taxonomy, gtdb_taxonomy)
                    elif gtdb_domain not in ncbi_taxonomy:
                        print '{0}\tgtdb_domain:{1}\tncbi_taxonomy:{2}\tgtdb_taxonomy:{3}'.format(genome_id, gtdb_domain, ncbi_taxonomy, gtdb_taxonomy)

            print "#End"

            # Compare NCBI and GTDB taxonomy
            query = ("SELECT g.id_at_source,mt.ncbi_taxonomy,gtv.gtdb_taxonomy FROM metadata_taxonomy mt " +
                     "LEFT JOIN genomes g using (id) LEFT JOIN gtdb_taxonomy_view gtv USING (id) " +
                     " WHERE genome_source_id in (2,3) ")
            self.cur.execute(query)
            list_domain = [[a, b, d] for (a, b, d) in self.cur]
            print "#Conflicting Taxonomy:"
            for genome_id, ncbi_taxonomy, gtdb_taxonomy in list_domain:
                if ncbi_taxonomy is None or gtdb_taxonomy is None:
                    continue

                for r in xrange(0, rank_depth + 1):
                    ncbi_taxon = ncbi_taxonomy.split(';')[r]
                    gtdb_taxon = gtdb_taxonomy.split(';')[r]
                    if ncbi_taxon == Taxonomy.rank_prefixes[r] or gtdb_taxon == Taxonomy.rank_prefixes[r]:
                        continue

                    if ncbi_taxon != gtdb_taxon:
                        print "{0}\tRank:{1}\tncbi_taxonomy:{2}\tgtdb_taxonomy:{3}".format(genome_id,
                                                                                           Taxonomy.rank_labels[r],
                                                                                           ncbi_taxonomy,
                                                                                           gtdb_taxonomy)
                        break
            print "#End"

        except GenomeDatabaseError as e:
            raise e

        return True

    def _chompRecord(self, record):
        if record.startswith("RS") or record.startswith("GB"):
            return record[3:]
        if record.startswith("U_"):
            return record[2:]

    def CheckUserIDsDuplicates(self):
        list_genome_ids = {}
        list_duplicates = []
        user_genome_dir = Config.GTDB_GENOME_USR_DIR
        for user_id in os.listdir(user_genome_dir):
            full_user_dir = os.path.join(user_genome_dir, user_id)
            if not os.path.isdir(full_user_dir):
                continue
            for genome_id in os.listdir(full_user_dir):
                if genome_id in list_genome_ids:
                    list_genome_ids[genome_id].append(full_user_dir)
                    list_duplicates.append(genome_id)
                else:
                    list_genome_ids[genome_id] = [full_user_dir]

        set_duplicates = set(list_duplicates)
        for dup in set_duplicates:
            print "Genome {0}:".format(dup)
            for path in list_genome_ids[dup]:
                print "- {0}".format(path)
        print "#############"

    def RunDomainConsistency(self):
        try:
            query = ("SELECT * from (SELECT id_at_source, " +
                     "CASE WHEN bac_mark::float/120*100 < 10 and arc_mark::float/122*100 < 10 THEN 'None' " +
                     "WHEN bac_mark::float/120*100 < arc_mark::float/122*100 THEN 'd__Archaea' " +
                     "ELSE 'd__Bacteria' END as marker_domain, " +
                     "gtdb_domain FROM genomes " +
                     "LEFT JOIN (SELECT id,count(*) as bac_mark from genomes g " +
                     "LEFT JOIN metadata_taxonomy USING (id) " +
                     "LEFT JOIN aligned_markers am on am.genome_id = g.id " +
                     "WHERE gtdb_representative is TRUE AND am.marker_id in (SELECT marker_id from marker_set_contents WHERE set_id =1) and am.evalue is not NULL " +
                     "group BY id) as tmpbac USING (id) " +
                     "LEFT JOIN (SELECT id,count(*) as arc_mark from genomes g " +
                     "LEFT JOIN metadata_taxonomy USING (id) " +
                     "LEFT JOIN aligned_markers am on am.genome_id = g.id " +
                     "WHERE gtdb_representative is TRUE AND am.marker_id in (SELECT marker_id from marker_set_contents WHERE set_id =2) and am.evalue is not NULL " +
                     "group BY id) as tmparc USING (id) " +
                     "LEFT JOIN metadata_taxonomy USING (id) " +
                     "WHERE gtdb_representative is TRUE ) as gtdb_difference " +
                     "WHERE marker_domain not like gtdb_domain")
            self.cur.execute(query)
            list_genome = [[a, b, d] for (a, b, d) in self.cur if b != d]
            if len(list_genome) > 0:
                print "Genome\tmarker_domain\tgtdb_domain"
                for item in list_genome:
                    print "{0}\t{1}\t{2}".format(item[0], item[1], item[2])
            print "Finished"
        except GenomeDatabaseError as e:
            raise e

        return True

    def RealignNCBIgenomes(self):
        try:
            query = ("SELECT g.id,COALESCE(marker_count,0) from genomes g " +
                     "LEFT JOIN (SELECT id_at_source,count(*) as marker_count from genomes g " +
                     "LEFT JOIN aligned_markers am ON am.genome_id = g.id " +
                     "LEFT JOIN marker_set_contents msc ON msc.marker_id = am.marker_id " +
                     "WHERE genome_source_id != 1 " +
                     "AND msc.set_id in (1,2) " +
                     "group by id_at_source) as marktmp USING (id_at_source) " +
                     "LEFT JOIN metadata_taxonomy mt on g.id = mt.id " +
                     "where g.genome_source_id != 1 and marker_count is NULL " +
                     "and gtdb_representative is not NULL " +
                     "ORDER BY marker_count")
            self.cur.execute(query)
            list_genome = [a for (a, _b) in self.cur]
            if len(list_genome) > 0:
                # get canonical bacterial and archaeal markers
                marker_set_mngr = MarkerSetManager(self.cur, self.currentUser)
                bac_marker_ids = marker_set_mngr.canonicalBacterialMarkers()
                ar_marker_ids = marker_set_mngr.canonicalArchaealMarkers()

                # identify and align genes from canonical bacterial and
                # archaeal marker sets
                all_markers = set(bac_marker_ids).union(ar_marker_ids)
                aligned_mngr = AlignedMarkerManager(
                    self.cur, self.threads, self.db_release)
                aligned_mngr.calculateAlignedMarkerSets(
                    list_genome, all_markers)
        except GenomeDatabaseError as e:
            raise e

        return True
