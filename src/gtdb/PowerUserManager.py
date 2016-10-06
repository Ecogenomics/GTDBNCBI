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


class PowerUserManager(object):

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

    def runTreeWeightedExceptions(self, path, comp, conta, qweight, qt):
        '''
        Function: runTreeWeightedException
        Export list of NCBI genomes that do comply the filter criteria but are of interest.

        :param path: Path to the output file
        '''
        try:
            if (not self.currentUser.isRootUser()):
                raise GenomeDatabaseError("Only the root user can run this command")
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
                            dict_except_order[ncbitax.split(';')[3]] = {'quality': float(qual), 'id': gid, 'full_info': [name, compl, conta, ncbitax, size, qual, orga, 'First']}
                        else:
                            if dict_except_order.get(ncbitax.split(';')[3]).get('quality') < float(qual):
                                dict_except_order[ncbitax.split(';')[3]] = {'quality': float(qual), 'id': gid, 'full_info': [name, compl, conta, ncbitax, size, qual, orga, 'First']}
                else:
                    dict_except_order[gid] = {'quality': float(qual), 'id': gid, 'full_info': [name, compl, conta, ncbitax, size, qual, orga, 'unknown order']}

            fh = open(path, "w")
            fh.write("Name,CheckM_Completeness,CheckM_Contamination,NCBI_Taxonomy,Genome_size,Quality_Threshold,Organism_name,Filter_passed\n")
            for _k, item in dict_except_order.iteritems():
                fh.write(",".join(str(v) for v in item.get('full_info')) + "\n")
            fh.close()

        except GenomeDatabaseError as e:
            raise e
        return True

    def _checkTaxonomyUniqueness(self, ncbitax, order_list):
        ncbi_tax_order = ncbitax.split(';')[3]
        if ncbi_tax_order in order_list:
            return False

        return True

    def runSanityCheck(self):
        try:
            if (not self.currentUser.isRootUser()):
                raise GenomeDatabaseError("Only the root user can run this command")
            query = ("SELECT id FROM genomes where genome_source_id in (2,3)")
            self.cur.execute(query)
            ncbi_ids = [gid for (gid,) in self.cur]

            query = ("SELECT id,id_at_source FROM genomes")
            self.cur.execute(query)
            raw_ids = [(gid, source_id) for (gid, source_id) in self.cur]
            all_ids, all_source_ids = zip(*raw_ids)
            dict_all_ids = {k: v for (k, v) in raw_ids}

            query = ("SELECT distinct(gtdb_genome_representative) from metadata_taxonomy where gtdb_genome_representative is not NULL")
            self.cur.execute(query)
            representatives = [self._chompRecord(record) for (record,) in self.cur]

            # check if the representatives are still in the database
            for representative in representatives:
                if representative not in all_source_ids:
                    print "REPRESENTATIVE {0} has been removed from the database".format(representative)

            query = ("SELECT id,ncbi_protein_count,ncbi_submitter from metadata_ncbi WHERE id in (SELECT id from genomes where genome_source_id in (2,3));")
            self.cur.execute(query)
            dict_meta_ncbi = {gid: {"count": count, "submitter": project} for (gid, count, project) in self.cur}

            for ncbi_genome in ncbi_ids:
                if ncbi_genome not in dict_meta_ncbi:
                    print "{0} has no metadata in metadata_ncbi".format(dict_all_ids[ncbi_genome])
                else:
                    if dict_meta_ncbi[ncbi_genome]["count"] is None or dict_meta_ncbi[ncbi_genome]["count"] == '' or dict_meta_ncbi[ncbi_genome]["count"] == 0:
                        print "{0} ncbi_protein_count value in metadata_ncbi is {1}".format(dict_all_ids[ncbi_genome], dict_meta_ncbi[ncbi_genome]["count"])
                    if dict_meta_ncbi[ncbi_genome]["submitter"] is None or dict_meta_ncbi[ncbi_genome]["submitter"] == '':
                        print "{0} ncbi_submitter value in metadata_ncbi is {1}".format(dict_all_ids[ncbi_genome], dict_meta_ncbi[ncbi_genome]["submitter"])

            query = ("SELECT id,checkm_completeness,protein_count from metadata_genes")
            self.cur.execute(query)
            dict_meta_genes = {gid: {"checkm": checkm, "protein_count": count} for (gid, checkm, count) in self.cur}

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
            # Check if gtdb_domain is the same as the gtdb_taxonomy and ncbi_taxonomy domain
            query = ("SELECT g.id_at_source, mt.gtdb_domain, mt.ncbi_taxonomy, mt.gtdb_taxonomy FROM metadata_taxonomy mt LEFT JOIN genomes g using (id) where genome_source_id in (2,3) ")
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
            query = ("SELECT g.id_at_source,mt.ncbi_taxonomy,mt.gtdb_taxonomy FROM metadata_taxonomy mt LEFT JOIN genomes g using (id) where genome_source_id in (2,3) ")
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
