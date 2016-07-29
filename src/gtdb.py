#!/usr/bin/env python

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

import argparse
import sys
import os
import inspect
import ntpath
import logging

from biolib.common import make_sure_path_exists
from biolib.external.execute import check_dependencies
from biolib.misc.custom_help_formatter import CustomHelpFormatter

from gtdb import GenomeDatabase
from gtdb import DefaultValues
from gtdb.Exceptions import GenomeDatabaseError


def version():
    """Read software and NCBI version information from file.

    Returns
    -------
    str
        Software version.
    str
        Software, NCBI database, and GTDB database versions.
    """
    cur_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    # cur_dir = os.path.dirname(os.path.realpath(__file__))
    version_file = open(os.path.join(cur_dir, 'VERSION'))

    software_version = version_file.readline().strip()
    software_version = software_version[software_version.find('=') + 1:]

    ncbi_version = version_file.readline().strip()
    ncbi_version = ncbi_version[ncbi_version.find('=') + 1:]

    gtdb_version = version_file.readline().strip()
    gtdb_version = gtdb_version[gtdb_version.find('=') + 1:]

    return software_version, ncbi_version, gtdb_version


def versionInfo():
    """Get version information.

    Returns
    -------
    str
        String indication software and NCBI version information.
    """

    software_version, ncbi_version, gtdb_version = version()
    return 'GTDB v%s (NCBI RefSeq %s; Internal database v%s)' % (software_version, ncbi_version, gtdb_version)


def loggerSetup(output_dir, silent=False):
    """Set logging for application.

    Parameters
    ----------
    output_dir : str
        Output directory for log file.
    silent : boolean
        Flag indicating if output to stdout should be suppressed.
    """

    # setup general properties of logger
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                   datefmt="%Y-%m-%d %H:%M:%S")

    # setup logging to console
    if not silent:
        stream_logger = logging.StreamHandler(sys.stdout)
        stream_logger.setFormatter(log_format)
        stream_logger.setLevel(logging.DEBUG)
        logger.addHandler(stream_logger)

    if output_dir:
        make_sure_path_exists(output_dir)
        file_logger = logging.FileHandler(
            os.path.join(output_dir, 'gtdb.log'), 'a')
        file_logger.setFormatter(log_format)
        logger.addHandler(file_logger)

    logger.info(versionInfo())
    logger.info(ntpath.basename(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]))


def DumpDBErrors(db):
    ErrorReport("\n".join(["\t" + x for x in db.GetErrors()]) + "\n")
    db.ClearErrors()


def DumpDBWarnings(db):
    ErrorReport("\n".join(["\t" + x for x in db.GetWarnings()]) + "\n")
    db.ClearWarnings()


def ErrorReport(msg):
    sys.stderr.write(msg)
    sys.stderr.flush()


def AddUser(db, args):
    log_has_root = False
    if args.login_as_root and args.has_root:
        log_has_root = True
    return db.addUser(args.username,
                      args.role,
                      log_has_root)


def EditUser(db, args):
    log_has_root = False
    if args.login_as_root and args.has_root:
        log_has_root = True
    return db.editUser(args.username,
                       args.role,
                       log_has_root)


def AddGenomes(db, args):

    return db.AddGenomes(args.batchfile,
                         args.checkm_file,
                         args.study_file,
                         args.genome_list_id,
                         args.genome_list_name)


def CreateTreeData(db, args):
    guaranteed_genomes = set()

    gid_list_and_guaranteed_gid = db.GetGenomeIds(args.all_dereplicated,
                                                  args.ncbi_dereplicated,
                                                  args.user_dereplicated,
                                                  args.donovan_sra_representatives,
                                                  args.all_genomes,
                                                  args.ncbi_genomes,
                                                  args.user_genomes,
                                                  args.genome_list_ids,
                                                  args.genome_ids,
                                                  args.genome_batchfile)
    genome_id_list = gid_list_and_guaranteed_gid[0]
    guaranteed_genomes.update(gid_list_and_guaranteed_gid[1])

    marker_id_list = db.GetMarkerIds(args.marker_ids,
                                     args.marker_set_ids,
                                     args.marker_batchfile)

    return db.MakeTreeData(marker_id_list, genome_id_list,
                           args.out_dir, args.prefix,
                           args.quality_threshold,
                           args.quality_weight,
                           args.comp_threshold,
                           args.cont_threshold,
                           args.min_perc_aa, args.min_perc_taxa, args.consensus,
                           args.taxa_filter,
                           args.excluded_genome_list_ids, args.excluded_genome_ids,
                           args.guaranteed_genome_list_ids,
                           args.guaranteed_genome_ids,
                           args.guaranteed_batchfile,
                           guaranteed_genomes,
                           not args.no_alignment,
                           args.individual,
                           not args.no_tree)


def ViewGenomes(db, args):
    if args.view_all:
        return db.ViewGenomes()
    else:
        external_ids = None
        if args.id_list:
            external_ids = args.id_list.split(",")
        return db.ViewGenomes(args.batchfile, external_ids)


def DeleteGenomes(db, args):

    external_ids = None
    if args.id_list:
        external_ids = args.id_list.split(",")

    list_ids = None
    if args.list_of_list_id:
        list_ids = args.list_of_list_id.split(",")

    return db.DeleteGenomes(args.batchfile, external_ids, list_ids, args.reason)


def PullGenomes(db, args):

    make_sure_path_exists(args.out_dir)

    external_ids = None
    if args.id_list:
        external_ids = args.id_list.split(",")

    list_ids = None
    if args.list_of_list_id:
        list_ids = args.list_of_list_id.split(",")

    return db.PullGenomes(args.batchfile,
                          external_ids,
                          list_ids,
                          args.genomic,
                          args.gene,
                          args.out_dir,
                          args.gtdb_header)


def ExportSSUSequences(db, args):
    return db.ExportSSUSequences(args.outfile)


def CreateGenomeList(db, args):

    external_ids = []

    if args.genome_ids:
        external_ids = args.genome_ids.split(",")
    genome_list_id = db.CreateGenomeList(args.batchfile,
                                         external_ids,
                                         args.name,
                                         args.description,
                                         (not args.public))

    if not genome_list_id:
        return False

    return genome_list_id


def ViewGenomeLists(db, args):

    genome_lists = []
    if args.root_owned or (args.self_owned and db.currentUser.isRootUser()):
        genome_lists = db.GetVisibleGenomeListsByOwner(include_private=False)
    elif args.self_owned:
        genome_lists = db.GetVisibleGenomeListsByOwner(
            owner_id=db.currentUser.getUserId(), include_private=False)
    elif args.all_public:
        genome_lists = db.GetAllVisibleGenomeListIds(include_private=False)
    elif args.all:
        genome_lists = db.GetAllVisibleGenomeListIds(include_private=True)
    elif args.owner_name is not None:
        genome_lists = db.GetGenomeListIdsforUser(args.owner_name)
    else:
        # TODO: this
        db.ReportError(
            "There is an unknown argument")
        return False

    if len(genome_lists) == 0:
        print "No genomes lists found."
        return True

    return db.PrintGenomeListsDetails(genome_lists)


def ContentsGenomeLists(db, args):

    list_ids = []
    if args.list_ids:
        list_ids = args.list_ids.split(",")

    return db.ViewGenomeListsContents(list_ids)


def EditGenomeLists(db, args):

    genome_ids = None
    if args.genome_ids:
        genome_ids = args.genome_ids.split(",")

    private = None
    if args.public:
        private = False
    if args.private:
        private = True

    return db.EditGenomeList(args.list_id, args.batchfile, genome_ids, args.operation, args.name, args.description, private)


def DeleteGenomeLists(db, args):
    list_ids = None
    if args.list_ids:
        list_ids = args.list_ids.split(",")
    return db.DeleteGenomeLists(list_ids)


def ViewMarkers(db, args):

    if args.view_all:
        return db.ViewMarkers()
    else:
        external_ids = None
        if args.id_list:
            external_ids = args.id_list.split(",")
        return db.ViewMarkers(args.batchfile, external_ids)


def CreateMarkerSet(db, args):

    external_marker_ids = []
    if args.marker_ids:
        external_marker_ids = args.marker_ids.split(",")

    marker_set_id = db.CreateMarkerSet(args.batchfile,
                                       external_marker_ids,
                                       args.name,
                                       args.description,
                                       (not args.public))

    return marker_set_id


def ViewMarkerSets(db, args):

    marker_sets = []
    if args.root_owned or (args.self_owned and db.currentUser.isRootUser()):
        marker_sets = db.GetVisibleMarkerSetsByOwner(include_private=False)
    elif args.self_owned:
        marker_sets = db.GetVisibleMarkerSetsByOwner(
            db.currentUser.getUserId(), include_private=False)
    elif args.all_public:
        marker_sets = db.GetAllVisibleMarkerSetIds(include_private=False)
    elif args.all:
        marker_sets = db.GetAllVisibleMarkerSetIds(include_private=True)
    elif args.owner_name is not None:
        marker_sets = db.GetAllMarkerSetsforUser(args.owner_name)
    else:
        # TODO: this
        db.ReportError(
            "There is an unknown argument")
        return False

    if len(marker_sets) == 0:
        print "No marker sets found."
        return True

    return db.printMarkerSetsDetails(marker_sets)


def EditMarkerSet(db, args):

    marker_ids = None
    if args.marker_ids:
        marker_ids = args.marker_ids.split(",")

    private = None
    if args.public:
        private = False
    if args.private:
        private = True

    return db.EditMarkerSet(args.set_id, args.batchfile, marker_ids, args.operation, args.name, args.description, private)


def DeleteMarkerSets(db, args):
    set_ids = None
    if args.set_ids:
        set_ids = args.set_ids.split(",")
    return db.deleteMarkerSets(set_ids)


def MarkerSetsContents(db, args):

    set_ids = []
    if args.set_ids:
        set_ids = args.set_ids.split(",")

    return db.ViewMarkerSetsContents(set_ids)


def viewMetadata(db, args):
    return db.ViewMetadata()


def exportMetadata(db, args):
    return db.ExportMetadata(args.outfile)


def importMetadata(db, args):
    if not db.currentUser.isRootUser():
        logging.getLogger().warning("Only the root user may import metadata.")
        return True

    return db.ImportMetadata(args.table, args.field, args.typemeta, args.metadatafile)


def createMetadata(db, args):
    if not db.currentUser.isRootUser():
        logging.getLogger().warning(
            "Only the root user may create new metadata fields.")

        return True

    return db.CreateMetadata(args.metadatafile)


def exportTaxonomyGTDB(db, args):
    return db.ExportTaxonomy('GTDB', args.outfile)


def exportTaxonomyNCBI(db, args):
    return db.ExportTaxonomy('NCBI', args.outfile)


def DatabaseStatsData(db, args):
    return db.ReportStats()


def RunTreeExceptions(db, args):
    return db.RunTreeWeightedExceptions(args.outfile,
                                        args.comp_threshold,
                                        args.cont_threshold,
                                        args.quality_weight,
                                        args.quality_threshold)


def RunSanityCheck(db, args):
    return db.RunSanityCheck()

def RunTaxonomyCheck(db, args):
    return db.RunTaxonomyCheck(args.rank_depth)

def RunDomainAssignmentReport(db, args):
    return db.RunDomainAssignmentReport(args.outfile)
    
def ExportGenomePaths(db, args):
    return db.ExportGenomePaths(args.outfile)


if __name__ == '__main__':
    # make sure all required dependencies are on the system path
    check_dependencies(
        ['prodigal', 'genometk', 'blastn', 'hmmsearch', 'pfam_search.pl'])

    # create the top-level parser
    parser = argparse.ArgumentParser(prog='gtdb', formatter_class=CustomHelpFormatter)
    parser.add_argument('-r', dest='login_as_root', action='store_true',
                        help='Login as the root user.')
    parser.add_argument('-u', dest='logon_as_user',
                        help='Logon as this user (implies -r).'),
    parser.add_argument('-t', dest='threads', type=int, default=1,
                        help='Maximum number of threads/cpus to use.')
    parser.add_argument('-f', dest='force', action='store_true',
                        help='Force all action (required to override warnings for certain actions).')
    parser.add_argument('-y', dest='assume_yes', action='store_true',
                        help='Assume yes to all confirmation prompts (useful for batch processing).')
    parser.add_argument('--tab_table', action='store_true',
                        help='Write tables as tab-separated.')
    parser.add_argument('--silent', action='store_true',
                        help='Suppress output to screen.')
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Run in debug mode.')
    parser.add_argument('--version', action='version', version=versionInfo(),
                        help='Show version information.')

    category_parser = parser.add_subparsers(
        dest='category_parser_name')

    user_category_parser = category_parser.add_parser('users',
                                                      formatter_class=CustomHelpFormatter,
                                                      help='Commands for adding and modifying users.')

    user_category_subparser = user_category_parser.add_subparsers(
        help='User command help.', dest='user_subparser_name')

    genome_category_parser = category_parser.add_parser('genomes',
                                                        formatter_class=CustomHelpFormatter,
                                                        help='Commands for adding, viewing, and removing genomes.')
    genome_category_subparser = genome_category_parser.add_subparsers(help='Genome command help.',
                                                                      dest='genome_subparser_name')

    genome_list_category_parser = category_parser.add_parser('genome_lists',
                                                             formatter_class=CustomHelpFormatter,
                                                             help='Commands for adding, viewing, and removing lists of genomes.')
    genome_list_category_subparser = genome_list_category_parser.add_subparsers(help='Genome list command help.',
                                                                                dest='genome_list_subparser_name')

    marker_category_parser = category_parser.add_parser('markers',
                                                        formatter_class=CustomHelpFormatter,
                                                        help='Commands for adding or viewing marker genes.')
    marker_category_subparser = marker_category_parser.add_subparsers(help='Marker command help.',
                                                                      dest='marker_subparser_name')

    marker_set_category_parser = category_parser.add_parser('marker_sets',
                                                            formatter_class=CustomHelpFormatter,
                                                            help='Commands for adding, viewing, and removing sets of markers.')
    marker_set_category_subparser = marker_set_category_parser.add_subparsers(help='Marker set command help.',
                                                                              dest='marker_sets_subparser_name')

    metadata_category_parser = category_parser.add_parser('metadata',
                                                          formatter_class=CustomHelpFormatter,
                                                          help='Commands for adding, viewing, and removing metadata fields and values.')
    metadata_category_subparser = metadata_category_parser.add_subparsers(help='Metadata command help.',
                                                                          dest='metadata_subparser_name')

    taxonomy_category_parser = category_parser.add_parser('taxonomy',
                                                          formatter_class=CustomHelpFormatter,
                                                          help='Commands for exporting the taxonomic information.')
    taxonomy_category_subparser = taxonomy_category_parser.add_subparsers(help='Taxonomy command help.',
                                                                          dest='taxonomy_subparser_name')

    tree_category_parser = category_parser.add_parser('tree',
                                                      formatter_class=CustomHelpFormatter,
                                                      help='Commands for inferring a phylogeny.')
    tree_category_subparser = tree_category_parser.add_subparsers(help='Tree command help.',
                                                                  dest='tree_subparser_name')

    db_stats_category_parser = category_parser.add_parser('db_stats',
                                                          formatter_class=CustomHelpFormatter,
                                                          help='Commands for viewing database statistics.')
    db_stats_category_subparser = db_stats_category_parser.add_subparsers(help='Database statistics command help.',
                                                                          dest='db_stats_subparser_name')

    power_category_parser = category_parser.add_parser('power',
                                                       formatter_class=CustomHelpFormatter,
                                                       help='Power command to run miscellaneous scripts')
    power_category_subparser = power_category_parser.add_subparsers(help='Power command help.',
                                                                    dest='power_subparser_name')

# -------- User Management subparsers

    # user add parser
    parser_user_add = user_category_subparser.add_parser('add',
                                                         add_help=False,
                                                         formatter_class=CustomHelpFormatter,
                                                         help='Add a user')

    required_user_add = parser_user_add.add_argument_group('required  arguments')
    required_user_add.add_argument('--username', dest='username', required=True,
                                   help='Username of the new user.')

    optional_user_add = parser_user_add.add_argument_group('optional arguments')
    optional_user_add.add_argument('--role', dest='role', choices=('user', 'admin'), required=False,
                                   help='Role of the new user.')
    optional_user_add.add_argument('--has_root', dest='has_root', action="store_true", required=False,
                                   help='User has permission to become the root user.')
    optional_user_add.add_argument('-h', '--help', action="help",
                                   help="Show help message.")

    parser_user_add.set_defaults(func=AddUser)

    # user edit parser
    parser_user_edit = user_category_subparser.add_parser('edit',
                                                          add_help=False,
                                                          formatter_class=CustomHelpFormatter,
                                                          help='Edit a user')

    required_user_edit = parser_user_edit.add_argument_group('required arguments')
    required_user_edit.add_argument('--username', dest='username', required=True,
                                    help='Username of the user to edit.')

    mutual_user_edit = parser_user_edit.add_argument_group('mutually exclusive optional arguments')
    mutex_group = mutual_user_edit.add_mutually_exclusive_group(required=False)
    mutex_group.add_argument('--has_root', dest='has_root', action="store_true", default=None,
                             help='Grant user the permission to become the root user.')
    mutex_group.add_argument('--no_root', dest='has_root', action="store_false", default=None,
                             help="Revoke user's permission to become the root user.")

    optional_user_edit = parser_user_edit.add_argument_group('optional arguments')
    optional_user_edit.add_argument('--role', dest='role', choices=('user', 'admin'), required=False,
                                    help='Change the user to this role.')
    optional_user_edit.add_argument('-h', '--help', action="help",
                                    help="Show help message.")

    parser_user_edit.set_defaults(func=EditUser)

# -------- Genome Management subparsers

    # genome add parser
    parser_genome_add = genome_category_subparser.add_parser('add',
                                                             add_help=False,
                                                             formatter_class=CustomHelpFormatter,
                                                             help='Add one or more genomes to the tree.')

    required_genome_add = parser_genome_add.add_argument_group('required arguments')
    required_genome_add.add_argument('--batchfile', dest='batchfile', required=True,
                                     help='Batch file describing genomes - one per line, tab separated in 3-6 columns (bin_filename, bin_name, bin_desc, [gene_filename], [source], [id_at_source]).')
    required_genome_add.add_argument('--checkm_results', dest='checkm_file', required=True,
                                     help='Provide a tab-separated CheckM results file (e.g. "checkm taxonomy_wf -f CHECKM_FILE --tab_table domain Bacteria ./bins output").')
    required_genome_add.add_argument('--study_file', required=True,
                                     help='File describing study and workflow from which genomes were recovered.')

    mutual_genome_add = parser_genome_add.add_argument_group('mutually exclusive required arguments')
    mutex_group = mutual_genome_add.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--modify_list', dest='genome_list_id',
                             help='Modify a genome list with the specified id and add all batchfile genomes into it.')
    mutex_group.add_argument('--create_list', dest='genome_list_name',
                             help='Create a genome list with the specified name and add all batchfile genomes into it.')
    mutex_group.add_argument('--no_list', dest='no_genome_list', action="store_true",
                             help="Don't add these genomes to a list.")

    optional_genome_add = parser_genome_add.add_argument_group('optional arguments')
    optional_genome_add.add_argument('-h', '--help', action="help",
                                     help="Show help message.")

    parser_genome_add.set_defaults(func=AddGenomes)

    # genome delete parser
    parser_genome_delete = genome_category_subparser.add_parser('delete',
                                                                add_help=False,
                                                                formatter_class=CustomHelpFormatter,
                                                                help='Remove genomes from the database.')
    atleastone_genome_delete = parser_genome_delete.add_argument_group('At least one argument required')
    atleastone_genome_delete.add_argument('--batchfile', dest='batchfile', default=None,
                                          help='Batchfile of genome IDs (one per line) to delete.')
    atleastone_genome_delete.add_argument('--genome_ids', dest='id_list', default=None,
                                          help='Provide a list of genome IDs (comma separated) to delete.')
    atleastone_genome_delete.add_argument('--list_ids', dest='list_of_list_id', default=None,
                                          help='Provide IDs of genome list (comma separated) to delete. Genomes part of those lists will be deleted and move to deprecated.')

    required_genome_delete = parser_genome_delete.add_argument_group('required named arguments')
    required_genome_delete.add_argument('--reason', dest='reason', required=True,
                                        help='Provide a reason why genomes are deleted.')

    optional_genome_delete = parser_genome_delete.add_argument_group('optional arguments')
    optional_genome_delete.add_argument('-h', '--help', action="help",
                                        help="Show help message.")

    parser_genome_delete.set_defaults(func=DeleteGenomes)

    # genome pull parser
    parser_genome_pull = genome_category_subparser.add_parser('pull',
                                                              add_help=False,
                                                              formatter_class=CustomHelpFormatter,
                                                              help='Pull genomic and gene data from database.')
    atleastone_genome_pull = parser_genome_pull.add_argument_group('At least one argument required')
    atleastone_genome_pull.add_argument('--batchfile', dest='batchfile', default=None,
                                        help='Batchfile of genome IDs (one per line) to pull.')
    atleastone_genome_pull.add_argument('--genome_ids', dest='id_list', default=None,
                                        help='Provide a list of genome IDs (comma separated) to pull.')
    atleastone_genome_pull.add_argument('--list_ids', dest='list_of_list_id', default=None,
                                        help='Provide IDs of genome lists (comma separated) to pull.')

    atleastone_genome_pull_data = parser_genome_pull.add_argument_group('At least one argument required')
    atleastone_genome_pull_data.add_argument('--genomic', default=None, action='store_true',
                                             help='Pull genomic sequences.')
    atleastone_genome_pull_data.add_argument('--gene', default=None, action='store_true',
                                             help='Pull called genes in amino acid space.')

    required_markers_genome_pull = parser_genome_pull.add_argument_group('required arguments')
    required_markers_genome_pull.add_argument('--output', dest='out_dir', required=True,
                                              help='Directory to output files.')

    optional_genome_pull = parser_genome_pull.add_argument_group('optional arguments')
    optional_genome_pull.add_argument('--gtdb_header', default=False, action="store_true",
                                      help="Add GTDB Prefix to NCBI Genomes (GB for Genbank and RS for Refseq).")
    optional_genome_pull.add_argument('-h', '--help', action="help",
                                      help="Show help message.")

    optional_genome_pull.set_defaults(func=PullGenomes)

    # genome view parser
    parser_genome_view = genome_category_subparser.add_parser('view',
                                                              add_help=False,
                                                              formatter_class=CustomHelpFormatter,
                                                              help='View the details of genomes in the database.')
    atleastone_genome_view = parser_genome_view.add_argument_group('At least one argument required')
    atleastone_genome_view.add_argument('--batchfile', dest='batchfile', default=None,
                                        help='Batchfile of genome IDs (one per line) to view.')
    atleastone_genome_view.add_argument('--genome_ids', dest='id_list', default=None,
                                        help='Provide a list of genome IDs (comma separated) to view.')
    atleastone_genome_view.add_argument('--all', dest='view_all', action="store_true",
                                        help='View all genomes in the database.')

    optional_genome_view = parser_genome_view.add_argument_group('optional arguments')
    optional_genome_view.add_argument('-h', '--help', action="help",
                                      help="Show help message.")

    parser_genome_view.set_defaults(func=ViewGenomes)

    # export ssu sequences for all genomes
    parser_genome_ssu_export = genome_category_subparser.add_parser('ssu_export',
                                                                    add_help=False,
                                                                    formatter_class=CustomHelpFormatter,
                                                                    help='Export a fasta file containing the SSU sequence best match for all genomes ')

    required_genome_ssu_export = parser_genome_ssu_export.add_argument_group('required arguments')
    required_genome_ssu_export.add_argument('--output', dest='outfile', default=None, required=True,
                                            help='Name of output file.')

    optional_genome_ssu_export = parser_genome_ssu_export.add_argument_group('optional arguments')
    optional_genome_ssu_export.add_argument('-h', '--help', action="help",
                                            help="Show help message.")

    parser_genome_ssu_export.set_defaults(func=ExportSSUSequences)

# -------- Genome Lists Management subparsers

    # Create genome list
    parser_gl_create = genome_list_category_subparser.add_parser('create',
                                                                 add_help=False,
                                                                 formatter_class=CustomHelpFormatter,
                                                                 help='Create a genome list.')
    required_gl_create = parser_gl_create.add_argument_group('required named arguments')
    required_gl_create.add_argument('--name', dest='name', required=True,
                                    help='Name of the genome list.')

    atleastone_gl_create = parser_gl_create.add_argument_group('At least one argument required')
    atleastone_gl_create.add_argument('--batchfile', dest='batchfile',
                                      help='File of genome IDs, one per line, to add to the create list.')
    atleastone_gl_create.add_argument('--genome_ids', dest='genome_ids',
                                      help='List of genome IDs (comma separated) to add to the create list.')

    optional_gl_create = parser_gl_create.add_argument_group('optional arguments')
    optional_gl_create.add_argument('--description', dest='description',
                                    help='A brief description of the genome list.')
    optional_gl_create.add_argument('--set_public', dest='public', action='store_true', default=False,
                                    help='Make the new list publicly visible.')
    optional_gl_create.add_argument('-h', '--help', action="help",
                                    help="Show help message.")

    parser_gl_create.set_defaults(func=CreateGenomeList)

    # View genome lists
    parser_gl_view = genome_list_category_subparser.add_parser('view',
                                                               add_help=False,
                                                               formatter_class=CustomHelpFormatter,
                                                               help='View genome lists.')

    mutual_genome_add = parser_gl_view.add_argument_group('mutually exclusive required arguments')
    mutex_group = mutual_genome_add.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--root', dest='root_owned', default=False, action='store_true',
                             help='Show genome lists owned by the root user.')
    mutex_group.add_argument('--self', dest='self_owned', default=False, action='store_true',
                             help='Show genome lists owned by you.')
    mutex_group.add_argument('--owner', dest='owner_name',
                             help='Show genome lists owned by a specific user.')
    mutex_group.add_argument('--all_public', default=False, action='store_true',
                             help='Show public genome lists from all users.')
    mutex_group.add_argument('--all', default=False, action='store_true',
                             help='View all genome lists.')

    optional_gl_view = parser_gl_view.add_argument_group('optional arguments')
    optional_gl_view.add_argument('-h', '--help', action="help",
                                  help="Show help message.")

    parser_gl_view.set_defaults(func=ViewGenomeLists)

    # Show genome list
    parser_gl_contents = genome_list_category_subparser.add_parser('contents',
                                                                   add_help=False,
                                                                   formatter_class=CustomHelpFormatter,
                                                                   help='View the contents of genome list(s).')
    required_gl_contents = parser_gl_contents.add_argument_group('required arguments')
    required_gl_contents.add_argument('--list_ids', dest='list_ids', required=True,
                                      help='Provide a list of genome list IDs (comma separated) whose contents you wish to view.')

    optional_gl_contents = parser_gl_contents.add_argument_group('optional arguments')
    optional_gl_contents.add_argument('-h', '--help', action="help",
                                      help="Show help message.")

    parser_gl_contents.set_defaults(func=ContentsGenomeLists)

    # Edit genome list
    parser_gl_edit = genome_list_category_subparser.add_parser('edit',
                                                               add_help=False,
                                                               formatter_class=CustomHelpFormatter,
                                                               help='Edit a genome list.')
    required_gl_edit = parser_gl_edit.add_argument_group('required named arguments')
    required_gl_edit.add_argument('--list_id', dest='list_id', required=True,
                                  help='Id of genome list to edit.')

    mutualoptional_gl_edit = parser_gl_edit.add_argument_group('mutually optional arguments')
    mutex_group = mutualoptional_gl_edit.add_mutually_exclusive_group(required=False)
    mutex_group.add_argument('--set_private', dest='private', action="store_true", default=False,
                             help='Make this genome list private (only you can see).')
    mutex_group.add_argument('--set_public', dest='public', action="store_true", default=False,
                             help='Make this genome list public (all users can see).')

    optional_gl_edit = parser_gl_edit.add_argument_group('optional arguments')
    optional_gl_edit.add_argument('--batchfile', dest='batchfile',
                                  help='A file of genome IDs, one per line, to add remove from the list.')
    optional_gl_edit.add_argument('--genome_ids', dest='genome_ids',
                                  help='List of genome IDs to add/remove from list.')
    optional_gl_edit.add_argument('--operation', dest='operation', choices=('add', 'remove'),
                                  help='What to do with the tree_ids with regards to the genome list. If all genomes are removed from a list, the list will be deleted.')
    optional_gl_edit.add_argument('--name', dest='name',
                                  help='Modify the name of the list to this.')
    optional_gl_edit.add_argument('--description', dest='description',
                                  help='Change the description of the genome list.')
    optional_gl_edit.add_argument('-h', '--help', action="help",
                                  help="Show help message.")

    parser_gl_edit.set_defaults(func=EditGenomeLists)

    # Delete genome list
    parser_gl_delete = genome_list_category_subparser.add_parser('delete',
                                                                 add_help=False,
                                                                 formatter_class=CustomHelpFormatter,
                                                                 help='Delete a genome list.')
    required_gl_delete = parser_gl_delete.add_argument_group('required named arguments')
    required_gl_delete.add_argument('--list_ids', dest='list_ids', required=True,
                                    help='Id of the genome lists to delete.')

    optional_gl_delete = parser_gl_delete.add_argument_group('optional arguments')
    optional_gl_delete.add_argument('-h', '--help', action="help",
                                    help="Show help message.")

    parser_gl_delete.set_defaults(func=DeleteGenomeLists)


# --------- Marker Management Subparsers

    # View markers
    parser_marker_view = marker_category_subparser.add_parser('view',
                                                              add_help=False,
                                                              formatter_class=CustomHelpFormatter,
                                                              help='View HMM markers in the database.')
    atleastone_marker_view = parser_marker_view.add_argument_group('At least one required argument')
    atleastone_marker_view.add_argument('--batchfile', dest='batchfile', default=None,
                                        help='Batchfile of marker IDs (one per line) to view.')
    atleastone_marker_view.add_argument('--marker_ids', dest='id_list', default=None,
                                        help='Provide a list of marker IDs (comma separated) to view.')
    atleastone_marker_view.add_argument('--all', dest='view_all', action="store_true",
                                        help='View all markers in the database.')

    optional_marker_view = parser_marker_view.add_argument_group('optional arguments')
    optional_marker_view.add_argument('-h', '--help', action="help",
                                      help="Show help message.")

    parser_marker_view.set_defaults(func=ViewMarkers)

# -------- Marker Set Management subparsers

    # Create marker set
    parser_ms_create = marker_set_category_subparser.add_parser('create',
                                                                add_help=False,
                                                                formatter_class=CustomHelpFormatter,
                                                                help='Create a marker set.')
    required_ms_create = parser_ms_create.add_argument_group('required arguments')
    required_ms_create.add_argument('--name', dest='name', required=True,
                                    help='Name of the marker set.')

    atleastone_ms_create = parser_ms_create.add_argument_group('At least one required argument')
    atleastone_ms_create.add_argument('--batchfile', dest='batchfile',
                                      help='File of marker IDs, one per line, to add to the created set.')
    atleastone_ms_create.add_argument('--marker_ids', dest='marker_ids',
                                      help='List of marker IDs (comma separated) to add to the created set.')

    optional_ms_create = parser_ms_create.add_argument_group('optional arguments')
    optional_ms_create.add_argument('--description', dest='description',
                                    help='Brief description of the marker set.')
    optional_ms_create.add_argument('--set_public', dest='public', action='store_true', default=False,
                                    help='Make the new set publicly visible.')
    optional_ms_create.add_argument('-h', '--help', action="help",
                                    help="Show help message.")

    parser_ms_create.set_defaults(func=CreateMarkerSet)

    parser_ms_view = marker_set_category_subparser.add_parser('view',
                                                              add_help=False,
                                                              formatter_class=CustomHelpFormatter,
                                                              help='View visible marker sets.')
    required_ms_view = parser_ms_view.add_argument_group('mutually exclusive required arguments')
    mutex_group = required_ms_view.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--root', dest='root_owned', default=False, action='store_true',
                             help='Only show marker sets owned by the root user.')
    mutex_group.add_argument('--self', dest='self_owned', default=False, action='store_true',
                             help='Only show marker sets owned by you.')
    mutex_group.add_argument('--owner', dest='owner_name',
                             help='Only show marker sets owned by a specific user.')
    mutex_group.add_argument('--all_public', default=False, action='store_true',
                             help='Show public marker sets from all users.')
    mutex_group.add_argument('--all', default=False, action='store_true',
                             help='View all marker sets.')

    optional_ms_view = parser_ms_view.add_argument_group('optional arguments')
    optional_ms_view.add_argument('-h', '--help', action="help",
                                  help="Show help message.")

    parser_ms_view.set_defaults(func=ViewMarkerSets)

    # Show marker set(s) contents
    parser_ms_contents = marker_set_category_subparser.add_parser('contents',
                                                                  add_help=False,
                                                                  formatter_class=CustomHelpFormatter,
                                                                  help='View the contents of marker set(s).')
    required_ms_contents = parser_ms_contents.add_argument_group('required named arguments')
    parser_ms_contents.add_argument('--set_ids', dest='set_ids', required=True,
                                    help='Provide a list of marker set IDs (comma separated) whose contents you wish to view.')

    optional_ms_contents = parser_ms_contents.add_argument_group('optional arguments')
    optional_ms_contents.add_argument('-h', '--help', action="help",
                                      help="Show help message.")

    parser_ms_contents.set_defaults(func=MarkerSetsContents)

    # Edit marker set
    parser_ms_edit = marker_set_category_subparser.add_parser('edit',
                                                              add_help=False,
                                                              formatter_class=CustomHelpFormatter,
                                                              help='Edit a marker set.')

    required_ms_edit = parser_ms_edit.add_argument_group('required arguments')
    required_ms_edit.add_argument('--set_id', dest='set_id', required=True,
                                  help='Id of the marker set to edit')

    mutual_ms_edit = parser_ms_edit.add_argument_group('mutually exclusive optional arguments')
    mutex_group = mutual_ms_edit.add_mutually_exclusive_group(required=False)
    mutex_group.add_argument('--set_private', dest='private', action="store_true", default=False,
                             help='Make this marker set private (only you can see).')
    mutex_group.add_argument('--set_public', dest='public', action="store_true", default=False,
                             help='Make this marker set public (all users can see).')

    optional_ms_edit = parser_ms_edit.add_argument_group('optional arguments')
    optional_ms_edit.add_argument('--batchfile', dest='batchfile',
                                  help='File of marker IDs, one per line, to add/remove from the set.')
    optional_ms_edit.add_argument('--marker_ids', dest='marker_ids',
                                  help='List (comma separated) of marker IDs to add/remove from set.')
    optional_ms_edit.add_argument('--operation', dest='operation', choices=('add', 'remove'),
                                  help='What to do with the provided marker IDs with regards to the marker set.')
    optional_ms_edit.add_argument('--name', dest='name',
                                  help='Modify the name of the set to this.')
    optional_ms_edit.add_argument('--description', dest='description',
                                  help='Change the description of the marker set to this.')
    optional_ms_edit.add_argument('-h', '--help', action="help",
                                  help="Show help message.")

    parser_ms_edit.set_defaults(func=EditMarkerSet)

    # Delete marker sets
    parser_ms_delete = marker_set_category_subparser.add_parser('delete',
                                                                add_help=False,
                                                                formatter_class=CustomHelpFormatter,
                                                                help='Delete a marker set.')
    required_ms_delete = parser_ms_delete.add_argument_group('required arguments')
    required_ms_delete.add_argument('--set_ids', dest='set_ids', required=True,
                                    help='List of marker set IDs (comma separated) whose contents you wish to delete.')

    optional_ms_delete = parser_ms_delete.add_argument_group('optional arguments')
    optional_ms_delete.add_argument('-h', '--help', action="help",
                                    help="Show help message.")

    parser_ms_delete.set_defaults(func=DeleteMarkerSets)

# -------- Metadata Management subparsers

    # metadata export parser
    parser_metadata_export = metadata_category_subparser.add_parser('export',
                                                                    add_help=False,
                                                                    formatter_class=CustomHelpFormatter,
                                                                    help='Export a CSV file with all metadata fields.')
    required_metadata_export = parser_metadata_export.add_argument_group('required arguments')
    required_metadata_export.add_argument('--output', dest='outfile', default=None, required=True,
                                          help='Name of output file.')

    optional_metadata_export = parser_metadata_export.add_argument_group('optional arguments')
    optional_metadata_export.add_argument('-h', '--help', action="help",
                                          help="Show help message.")

    parser_metadata_export.set_defaults(func=exportMetadata)

    # metadata view parser
    parser_metadata_view = metadata_category_subparser.add_parser('view',
                                                                  add_help=False,
                                                                  formatter_class=CustomHelpFormatter,
                                                                  help='List existing metadata fields with table name and description.')

    optional_metadata_view = parser_metadata_view.add_argument_group('optional arguments')
    optional_metadata_view.add_argument('-h', '--help', action="help",
                                        help="Show help message.")

    parser_metadata_view.set_defaults(func=viewMetadata)

    # metadata create columns parser
    parser_metadata_create = metadata_category_subparser.add_parser('create',
                                                                    add_help=False,
                                                                    formatter_class=CustomHelpFormatter,
                                                                    help='Create one or more new metadata field.')
    required_metadata_create = parser_metadata_create.add_argument_group('required arguments')
    required_metadata_create.add_argument('--file', dest='metadatafile',
                                          required=True, help='Metadata file describing the new fields - ' +
                                          'one field per line, tab separated in 4 columns' +
                                          '(name, description, datatype, metadata table).')

    optional_metadata_create = parser_metadata_create.add_argument_group('optional arguments')
    optional_metadata_create.add_argument('-h', '--help', action="help",
                                          help="Show help message.")

    parser_metadata_create.set_defaults(func=createMetadata)

    # metadata import parser
    parser_metadata_import = metadata_category_subparser.add_parser('import',
                                                                    add_help=False,
                                                                    formatter_class=CustomHelpFormatter,
                                                                    help='Import metadata values for a list of genome.')
    required_metadata_import = parser_metadata_import.add_argument_group('required arguments')
    required_metadata_import.add_argument('--table', dest='table', default=None, required=True,
                                          help='Table where the metadata field is present.')
    required_metadata_import.add_argument('--field', dest='field', default=None, required=True,
                                          help='Metadata field where the value(s) will be saved.')
    required_metadata_import.add_argument('--type', dest='typemeta', default=None, required=True,
                                          help='Type of the Metadata field.')
    required_metadata_import.add_argument('--metadatafile', dest='metadatafile', default=None, required=True,
                                          help='TSV file. One genome per line, tab separated in 2 columns indicating genome id and metadata value.')

    optional_metadata_import = parser_metadata_import.add_argument_group('optional arguments')
    optional_metadata_import.add_argument('-h', '--help', action="help",
                                          help="Show help message.")

    parser_metadata_import.set_defaults(func=importMetadata)

# -------- Taxonomy subparsers

    # GTDB taxonmoy export parser
    taxonomy_gtdb_export = taxonomy_category_subparser.add_parser('gtdb_export',
                                                                  add_help=False,
                                                                  formatter_class=CustomHelpFormatter,
                                                                  help='Export GTDB taxonomy as a TSV file.')
    required_taxonomy_gtdb_export = taxonomy_gtdb_export.add_argument_group('required arguments')
    required_taxonomy_gtdb_export.add_argument('--output', dest='outfile', default=None, required=True,
                                               help='Name of output file.')

    optional_taxonomy_gtdb_export = taxonomy_gtdb_export.add_argument_group('optional arguments')
    optional_taxonomy_gtdb_export.add_argument('-h', '--help', action="help",
                                               help="Show help message.")

    taxonomy_gtdb_export.set_defaults(func=exportTaxonomyGTDB)

    # NCBI taxonmoy export parser
    taxonomy_ncbi_export = taxonomy_category_subparser.add_parser('ncbi_export',
                                                                  add_help=False,
                                                                  formatter_class=CustomHelpFormatter,
                                                                  help='Export NCBI taxonomy as a TSV file.')
    required_taxonomy_ncbi_export = taxonomy_ncbi_export.add_argument_group('required arguments')
    required_taxonomy_ncbi_export.add_argument('--output', dest='outfile', default=None, required=True,
                                               help='Name of output file.')

    optional_taxonomy_ncbi_export = taxonomy_ncbi_export.add_argument_group('optional arguments')
    optional_taxonomy_ncbi_export.add_argument('-h', '--help', action="help",
                                               help="Show help message.")

    taxonomy_ncbi_export.set_defaults(func=exportTaxonomyNCBI)

# -------- Generate Tree Data
    parser_tree_create = tree_category_subparser.add_parser('create',
                                                            add_help=False,
                                                            formatter_class=CustomHelpFormatter,
                                                            help='Infer genome tree.')

    atleastone_genomes_create_tree = parser_tree_create.add_argument_group('minimum of one argument required')
    atleastone_genomes_create_tree.add_argument('--all_dereplicated', default=False, action='store_true',
                                                help=('Include all representative genomes and all genomes without ' +
                                                      'a representative. This is the set of genomes typically used ' +
                                                      'for initial tree exploration. Genomes are subject ' +
                                                      'to filtering.'))
    atleastone_genomes_create_tree.add_argument('--ncbi_dereplicated', default=False, action='store_true',
                                                help=('Include NCBI representative genomes and NCBI genomes without ' +
                                                      'a representative. This is the set of genomes typically used ' +
                                                      'for trees being published. Genomes are subject ' +
                                                      'to filtering.'))
    atleastone_genomes_create_tree.add_argument('--user_dereplicated', default=False, action='store_true',
                                                help=('Include User representative genomes and User genomes without ' +
                                                      'a representative. Genomes are subject to filtering.'))

    atleastone_genomes_create_tree.add_argument('--donovan_sra_representatives', default=False, action='store_true',
                                                help=('Include SRA representative genomes generated from Donovan SRA bins.' +
                                                      ' This is a temporary flag.'))

    atleastone_genomes_create_tree.add_argument('--all_genomes', default=False, action='store_true',
                                                help='Include all genomes, subject to filtering.')
    atleastone_genomes_create_tree.add_argument('--ncbi_genomes', default=False, action='store_true',
                                                help='Include all NCBI genomes, subject to filtering.')
    atleastone_genomes_create_tree.add_argument('--user_genomes', default=False, action='store_true',
                                                help='Include all User genomes, subject to filtering.')

    atleastone_genomes_create_tree.add_argument('--genome_list_ids', dest='genome_list_ids', default=None,
                                                help=('Genome list IDs (comma separated) whose ' +
                                                      'genomes should be included in the tree. Genomes are ' +
                                                      'subject to filtering.'))
    atleastone_genomes_create_tree.add_argument('--genome_ids', dest='genome_ids', default=None,
                                                help=('Genome IDs (comma separated) whose ' +
                                                      'genomes should be included in the tree. Genomes ' +
                                                      'are subject to filtering.'))
    atleastone_genomes_create_tree.add_argument('--genome_batchfile', dest='genome_batchfile', default=None,
                                                help=('File of genome IDs, one per line, to include in ' +
                                                      'the tree. Genomes are subject to filtering.'))

    atleastone_markers_create_tree = parser_tree_create.add_argument_group('minimum of one argument required')
    atleastone_markers_create_tree.add_argument('--marker_set_ids', dest='marker_set_ids', default=None,
                                                help='Marker set IDs (comma separated) whose markers will be used to infer the tree.')
    atleastone_markers_create_tree.add_argument('--marker_ids', dest='marker_ids', default=None,
                                                help='Marker IDs (comma separated) of markers that will be used to infer the tree.')
    atleastone_markers_create_tree.add_argument('--marker_batchfile', dest='marker_batchfile', default=None,
                                                help='File of marker IDs, one per line, to use when inferring the tree.')

    required_markers_create_tree = parser_tree_create.add_argument_group('required arguments')
    required_markers_create_tree.add_argument('--output', dest='out_dir', required=True,
                                              help='Directory to output files.')

    optional_markers_create_tree = parser_tree_create.add_argument_group('optional arguments')
    optional_markers_create_tree.add_argument('--quality_threshold', type=float, default=DefaultValues.DEFAULT_QUALITY_THRESHOLD,
                                              help='Filter genomes with a quality (completeness - weight*contamination) below threshold.')
    optional_markers_create_tree.add_argument('--quality_weight', type=float, default=DefaultValues.DEFAULT_QUALITY_WEIGHT,
                                              help='Weighting used to assess genome quality (completeness - weight*contamination).')

    optional_markers_create_tree.add_argument('--completeness_threshold', dest='comp_threshold', type=float, default=DefaultValues.DEFAULT_CHECKM_COMPLETENESS,
                                              help='Filter genomes below completeness threshold.')
    optional_markers_create_tree.add_argument('--contamination_threshold', dest='cont_threshold', type=float, default=DefaultValues.DEFAULT_CHECKM_CONTAMINATION,
                                              help='Filter genomes above contamination threshold.')

    optional_markers_create_tree.add_argument('--min_perc_aa', type=float, default=50,
                                              help='Filter genomes with an insufficient percentage of AA in the MSA.')
    optional_markers_create_tree.add_argument('--min_perc_taxa', type=float, default=50,
                                              help='minimum percentage of taxa required required to retain column.')
    optional_markers_create_tree.add_argument('--consensus', type=float, default=25,
                                              help='minimum percentage of the same amino acid required to retain column.')
    optional_markers_create_tree.add_argument('--excluded_genome_list_ids',
                                              help='Genome list IDs (comma separated) indicating genomes to exclude from the tree.')
    optional_markers_create_tree.add_argument('--excluded_genome_ids',
                                              help='Genome IDs (comma separated) indicating genomes to exclude from the tree.')
    optional_markers_create_tree.add_argument('--guaranteed_genome_list_ids',
                                              help='Genome list IDs (comma separated) indicating genomes to retain in the tree independent of filtering criteria.')
    optional_markers_create_tree.add_argument('--guaranteed_genome_ids',
                                              help='Genome IDs (comma separated) indicating genomes to retain in the tree independent of filtering criteria.')
    optional_markers_create_tree.add_argument('--guaranteed_batchfile',
                                              help='File of genome IDs, one per line, indicating genomes to retain in the tree independent of filtering criteria.')

    optional_markers_create_tree.add_argument('--taxa_filter',
                                              help='Filter genomes to taxa (comma separated) within specific taxonomic groups (e.g., d__Archaea or p__Proteobacteria, p__Actinobacteria).')

    optional_markers_create_tree.add_argument('--prefix', required=False, default='gtdb',
                                              help='Desired prefix for output files.')
    optional_markers_create_tree.add_argument('--no_alignment', action='store_true',
                                              help='Remove concatenated alignment in ARB metadata file.')
    optional_markers_create_tree.add_argument('--individual', action='store_true',
                                              help='Create individual FASTA files for each marker.')

    optional_markers_create_tree.add_argument('--no_tree', dest='no_tree', action="store_true",
                                              help="Output tree data, but do not infer a tree.")
    optional_markers_create_tree.add_argument('-h', '--help', action="help",
                                              help="Show help message.")

    parser_tree_create.set_defaults(func=CreateTreeData)

# -------- Generate Tree Data
    parser_db_stats_view = db_stats_category_subparser.add_parser('view',
                                                                  add_help=False,
                                                                  formatter_class=CustomHelpFormatter,
                                                                  help='View database statistics.')

    optional_db_stats_view = parser_db_stats_view.add_argument_group('optional arguments')
    optional_db_stats_view.add_argument('-h', '--help', action="help",
                                        help="Show help message.")

    parser_db_stats_view.set_defaults(func=DatabaseStatsData)


# Power user
# -------- Find exceptions in the tree
    parser_power_tree_exception = power_category_subparser.add_parser('weighted_tree_exceptions',
                                                                      add_help=False,
                                                                      formatter_class=CustomHelpFormatter,
                                                                      help='View NCBI records where genus is not present in the final tree and the checkm_completeness and contamination are less significant than the default values ')

    required_power_tree_exceptions = parser_power_tree_exception.add_argument_group('required arguments')
    required_power_tree_exceptions.add_argument('--output', dest='outfile', default=None, required=True,
                                                help='Name of output file.')

    optional_power_view = parser_power_tree_exception.add_argument_group('optional arguments')
    optional_power_view.add_argument('--quality_threshold', type=float, default=DefaultValues.DEFAULT_QUALITY_THRESHOLD,
                                     help='Filter genomes with a quality (completeness - weight*contamination) below threshold.')
    optional_power_view.add_argument('--quality_weight', type=float, default=DefaultValues.DEFAULT_QUALITY_WEIGHT,
                                     help='Weighting used to assess genome quality (completeness - weight*contamination).')
    optional_power_view.add_argument('--completeness_threshold', dest='comp_threshold', type=float, default=DefaultValues.DEFAULT_CHECKM_COMPLETENESS,
                                     help='Filter genomes below completeness threshold.')
    optional_power_view.add_argument('--contamination_threshold', dest='cont_threshold', type=float, default=DefaultValues.DEFAULT_CHECKM_CONTAMINATION,
                                     help='Filter genomes above contamination threshold.')
    optional_power_view.add_argument('-h', '--help', action="help",
                                     help="Show help message.")

    parser_power_tree_exception.set_defaults(func=RunTreeExceptions)

# --------- Export genome folder paths

    parser_power_genome_path = power_category_subparser.add_parser('genome_paths',
                                                                   add_help=False,
                                                                   formatter_class=CustomHelpFormatter,
                                                                   help='Export all genome path options')

    required_power_tree_genome_path = parser_power_genome_path.add_argument_group('required arguments')
    required_power_tree_genome_path.add_argument('--output', dest='outfile', default=None, required=True,
                                                 help='Name of output file.')

    optional_power_genome_path_export = parser_power_genome_path.add_argument_group('optional arguments')
    optional_power_genome_path_export.add_argument('-h', '--help', action="help",
                                                   help="Show help message.")

    parser_power_genome_path.set_defaults(func=ExportGenomePaths)

    # -------- Sanity check
    parser_sanity_exception = power_category_subparser.add_parser('sanity_check',
                                                                  add_help=False,
                                                                  formatter_class=CustomHelpFormatter,
                                                                  help='Run some sanity checks to see if all records are properly stored')

    optional_sanity_view = parser_sanity_exception.add_argument_group('optional arguments')
    optional_sanity_view.add_argument('-h', '--help', action="help",
                                      help="Show help message.")

    parser_sanity_exception.set_defaults(func=RunSanityCheck)

    # -------- Taxonomy check
    parser_taxonomy_check = power_category_subparser.add_parser('taxonomy_check',
                                                                  add_help=False,
                                                                  formatter_class=CustomHelpFormatter,
                                                                  help='Compare GTDB to NCBI taxonomy and report differences.')

    optional_taxonomy_check = parser_taxonomy_check.add_argument_group('optional arguments')
    optional_taxonomy_check.add_argument('--rank_depth', type=int, default=0,
                                      help='Deepest taxonomic rank to check: 0 (domain) to 6 (species).')
    optional_taxonomy_check.add_argument('-h', '--help', action="help",
                                      help="Show help message.")

    parser_taxonomy_check.set_defaults(func=RunTaxonomyCheck)
    
    # -------- Taxonomy check
    parser_domain_report = power_category_subparser.add_parser('domain_report',
                                                                  add_help=False,
                                                                  formatter_class=CustomHelpFormatter,
                                                                  help='Reports results of automated domain assignment.')

    required_domain_report = parser_domain_report.add_argument_group('required arguments')
    required_domain_report.add_argument('--output', dest='outfile', default=None, required=True,
                                                 help='Name of output file.')
                                                 
    optional_domain_report = parser_domain_report.add_argument_group('optional arguments')
    optional_domain_report.add_argument('-h', '--help', action="help",
                                      help="Show help message.")

    parser_domain_report.set_defaults(func=RunDomainAssignmentReport)

    # Do the parsing
    args = parser.parse_args()

    # setup logger
    if hasattr(args, 'out_dir'):
        loggerSetup(args.out_dir, args.silent)
    else:
        loggerSetup(None, args.silent)

    # Special parser checks
    if (args.category_parser_name == 'tree' and args.tree_subparser_name == 'create'):
        if (not args.all_dereplicated and
                not args.ncbi_dereplicated and
                not args.user_dereplicated and
                not args.donovan_sra_representatives and
                not args.all_genomes and
                not args.ncbi_genomes and
                not args.user_genomes and
                not args.genome_list_ids and
                not args.genome_ids and
                not args.genome_batchfile):
            parser_tree_create.error(
                'Need to specify at least one of --all_dereplicated, --ncbi_dereplicated, --user_dereplicated, --donovan_sra_representatives, --all_genomes, --ncbi_genomes, --user_genomes --genome_list_ids, --genome_ids, or --genome_batchfile.')

        if (not args.marker_set_ids
                and not args.marker_ids
                and not args.marker_batchfile):
            parser_tree_create.error(
                'Need to specify at least one of --marker_set_ids, --marker_ids or --marker_batchfile.')

    if (args.category_parser_name == 'genomes' and args.genome_subparser_name == 'view'):
        if (args.batchfile is not None or args.id_list is not None):
            if args.view_all:
                parser_genome_view.error(
                    'Argument --all must be used by itself.')
        elif not args.view_all:
            parser_genome_view.error(
                'Need to specify at least one of --all, --batchfile or --genome_ids.')

    if (args.category_parser_name == 'genomes' and args.genome_subparser_name == 'delete'):
        if (args.batchfile is None and args.id_list is None and args.list_of_list_id is None):
            parser_genome_delete.error(
                'Need to specify at least one of --batchfile or --genome_ids or --list_ids.')

    if (args.category_parser_name == 'markers' and args.marker_subparser_name == 'view'):
        if (args.batchfile is not None or args.id_list is not None):
            if args.view_all:
                parser_marker_view.error(
                    'Argument --all must be used by itself.')
        elif not args.view_all:
            parser_marker_view.error(
                'Need to specify at least one of --all, --batchfile or --marker_ids.')

    # initialise the backend
    db = GenomeDatabase.GenomeDatabase(args.threads, args.tab_table)
    db.conn.MakePostgresConnection()

    if args.debug:
        db.SetDebugMode(True)

    # Login
    try:
        db.Login(args.logon_as_user, args.login_as_root)
    except GenomeDatabaseError as e:
        db.conn.ClosePostgresConnection()
        ErrorReport(e.message + " The following error(s) were reported:\n")
        DumpDBErrors(db)
        sys.exit(-1)

    try:
        result = args.func(db, args)
    except:
        db.conn.ClosePostgresConnection()
        ErrorReport("Exception caught. Dumping info.\n")

        if db.GetWarnings():
            ErrorReport("Database reported the following warning(s):\n")
            DumpDBWarnings(db)

        if db.GetErrors():
            ErrorReport("Database reported the following errors(s):\n")
            DumpDBErrors(db)

        raise

    if db.GetWarnings():
        ErrorReport("Database reported the following warning(s):\n")
        DumpDBWarnings(db)

    if not result:
        ErrorReport(
            "Database action failed. The following error(s) were reported:\n")
        DumpDBErrors(db)

    db.conn.ClosePostgresConnection()
