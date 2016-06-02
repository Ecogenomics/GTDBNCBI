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
from Exceptions import GenomeDatabaseError


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

    def runTreeWeightedExceptions(self, path):
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
                             "WHERE g.genome_source_id IN (2,3) and " +
                             "mg.checkm_completeness > %s and mg.checkm_contamination < %s " +
                             "and mg.checkm_completeness-4*mg.checkm_contamination > %s and mt.ncbi_taxonomy is not NULL",
                             (DefaultValues.DEFAULT_CHECKM_COMPLETENESS, DefaultValues.DEFAULT_CHECKM_CONTAMINATION, DefaultValues.DEFAULT_QUALITY_THRESHOLD))
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

            dict_except_genus = {}
            for (gid, name, compl, conta, ncbitax, size, qual, orga) in self.cur:
                if self._checkTaxonomyUniqueness(ncbitax, order_list):
                    if self._checkTaxonomyUniqueness(ncbitax, dict_except_genus):
                        dict_except_genus[ncbitax] = {'quality': float(qual), 'id': gid, 'full_info': [name, compl, conta, ncbitax, size, qual, orga]}
                    else:
                        if dict_except_genus.get(ncbitax).get('quality') > float(qual):
                            dict_except_genus[ncbitax] = {'quality': float(qual), 'id': gid, 'full_info': [name, compl, conta, ncbitax, size, qual, orga]}

            list_exception_id = [info.get('id') for k, info in dict_except_genus.iteritems()]
            combined_list = list_exception_id + list(existing_id)
            self.cur.execute("SELECT g.name,mg.checkm_completeness,mg.checkm_contamination,mt.ncbi_taxonomy,mnuc.genome_size,(mg.checkm_completeness-4*mg.checkm_contamination) as quality_threshold,mn.ncbi_organism_name " +
                             "FROM genomes g " +
                             "LEFT JOIN metadata_genes mg USING (id) " +
                             "LEFT JOIN metadata_ncbi mn  USING (id) " +
                             "LEFT JOIN metadata_nucleotide mnuc  USING (id) " +
                             "LEFT JOIN metadata_taxonomy mt  USING (id) " +
                             "WHERE g.genome_source_id IN (2,3) and " +
                             "(mg.checkm_completeness > %s and  mg.checkm_contamination < %s " +
                             "and mg.checkm_completeness-4*mg.checkm_contamination > %s) and g.id not in %s",
                             (DefaultValues.EXCEPTION_FILTER_TWO_CHECKM_COMPLETENESS, DefaultValues.EXCEPTION_FILTER_TWO_CHECKM_CONTAMINATION, DefaultValues.EXCEPTION_FILTER_TWO_QUALITY_THRESHOLD, tuple(combined_list)))
            exception_2nd_filter = [[name, compl, conta, ncbitax, size, qual, orga] for (name, compl, conta, ncbitax, size, qual, orga) in self.cur]

            fh = open(path, "w")
            fh.write("Name,CheckM_Completeness,CheckM_Contamination,NCBI_Taxonomy,Genome_size,Quality_Threshold,Organism_name\n")
            for k, item in dict_except_genus.iteritems():
                fh.write(",".join(str(v) for v in item.get('full_info')) + "\n")
            fh.close
            for item in exception_2nd_filter:
                fh.write(",".join(str(v) for v in item) + "\n")
            fh.close

        except GenomeDatabaseError as e:
            raise e
        return True

    def _checkTaxonomyUniqueness(self, ncbitax, order_list):
        print '_checkTaxonomyUniqueness'
        ncbi_tax_order = ncbitax.split(';')[3]
        if ncbi_tax_order in order_list:
            return False

        return True

    def runTreeExceptions(self, path, filtered):
        '''
        Function: RunTreeException
        Export excluded NCBI records for a tree creation (with all default parameters) to a csv file

        :param path: Path to the output file
        '''
        try:
            if (not self.currentUser.isRootUser()):
                raise GenomeDatabaseError("Only the root user can run this command")
            self.cur.execute("SELECT distinct(string_to_array(mn.ncbi_organism_name, ' '))[1] " +
                             "FROM genomes g " +
                             "LEFT JOIN metadata_genes mg USING (id) " +
                             "LEFT JOIN metadata_ncbi mn  USING (id) " +
                             "WHERE g.genome_source_id IN (2,3) and " +
                             "mg.checkm_completeness > %s and mg.checkm_contamination < %s " +
                             "and mg.checkm_completeness-4*mg.checkm_contamination > %s and mn.ncbi_organism_name is not NULL",
                             (DefaultValues.DEFAULT_CHECKM_COMPLETENESS, DefaultValues.DEFAULT_CHECKM_CONTAMINATION, DefaultValues.DEFAULT_QUALITY_THRESHOLD))
            existing_genus = [genus for (genus,) in self.cur]
            self.cur.execute("SELECT g.name,mg.checkm_completeness,mg.checkm_contamination,mt.ncbi_taxonomy,mnuc.genome_size,(mg.checkm_completeness-4*mg.checkm_contamination) as quality_threshold,mn.ncbi_organism_name " +
                             "FROM genomes g " +
                             "LEFT JOIN metadata_genes mg USING (id) " +
                             "LEFT JOIN metadata_ncbi mn  USING (id) " +
                             "LEFT JOIN metadata_nucleotide mnuc  USING (id) " +
                             "LEFT JOIN metadata_taxonomy mt  USING (id) " +
                             "WHERE g.genome_source_id IN (2,3) and " +
                             "(mg.checkm_completeness < %s or mg.checkm_contamination > %s " +
                             "or mg.checkm_completeness-4*mg.checkm_contamination < %s) and mn.ncbi_organism_name is not NULL",
                             (DefaultValues.DEFAULT_CHECKM_COMPLETENESS, DefaultValues.DEFAULT_CHECKM_CONTAMINATION, DefaultValues.DEFAULT_QUALITY_THRESHOLD))

            fh = open(path, "w")
            fh.write("Name,CheckM_Completeness,CheckM_Contamination,NCBI_Taxonomy,Genome_size,Quality_Threshold,Organism_name\n")
            if filtered:
                exception_genus = [[name, compl, conta, ncbitax, size, qual, self._checkCandidatusPresence(orga)] for (name, compl, conta, ncbitax, size, qual, orga) in self.cur if self._filterOrganismName(orga, existing_genus)]
            else:
                exception_genus = [[name, compl, conta, ncbitax, size, qual, orga] for (name, compl, conta, ncbitax, size, qual, orga) in self.cur]
            for item in exception_genus:
                fh.write(",".join(str(v) for v in item) + "\n")
            fh.close

        except GenomeDatabaseError as e:
            raise e
        return True

    def runSanityCheck(self, path_file=None):
        try:
            if (not self.currentUser.isRootUser()):
                raise GenomeDatabaseError("Only the root user can run this command")
            query = ("SELECT id FROM genomes where genome_source_id in (2,3)")
            self.cur.execute(query)
            ncbi_ids = [id for (id,) in self.cur]

            query = ("SELECT id,id_at_source FROM genomes")
            self.cur.execute(query)
            raw_ids = [(id, source_id) for (id, source_id) in self.cur]
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
            dict_meta_ncbi = {id: {"count": count, "submitter": project} for (id, count, project) in self.cur}

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
            dict_meta_genes = {id: {"checkm": checkm, "protein_count": count} for (id, checkm, count) in self.cur}

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
            dict_meta_nuc = {id: {"gc": gc} for (id, gc) in self.cur}

            for genome in all_ids:
                if genome not in dict_meta_nuc:
                    print "{0} has no metadata in metadata_nucleotide".format(dict_all_ids[genome])
                else:
                    if dict_meta_nuc[genome]["gc"] is None or dict_meta_nuc[genome]["gc"] == '' or dict_meta_nuc[genome]["gc"] == 0:
                        print "{0} gc_count value in metadata_nucleotide is {1}".format(dict_all_ids[genome], dict_meta_nuc[genome]["gc"])
        except GenomeDatabaseError as e:
            raise e
        return True

    def _chompRecord(self, record):
        if record.startswith("RS") or record.startswith("GB"):
            return record[3:]
        if record.startswith("U_"):
            return record[2:]

    def _filterOrganismName(self, orga_name, existing_genus):
        if " BACTERIUM " in orga_name.upper():
            return False
        orga_name = self._checkCandidatusPresence(orga_name)
        if " sp. " in orga_name:
            return False
        if orga_name[0].islower():
            return False
        if len(orga_name.split(" ")) < 2:
            return False
        elif orga_name.split(" ")[0] in existing_genus:
            return False
        return True

    def _checkCandidatusPresence(self, orga_name):
        if "CANDIDATUS" in orga_name.upper():
            full_name_list = orga_name.split(" ")
            for element in full_name_list:
                if "CANDIDATUS" == element.upper():
                    full_name_list.remove(element)
            orga_name = " ".join(full_name_list)
        return orga_name
