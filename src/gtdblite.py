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
import pwd

from gtdblite import GenomeDatabase
from gtdblite import profiles
from gtdblite.Exceptions import GenomeDatabaseError


def GetLinuxUsername():
    return pwd.getpwuid(os.getuid())[0]


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
    has_root = False
    if args.has_root:
        has_root = True

    return db.AddUser(args.username, args.role, has_root)


def EditUser(db, args):
    return db.EditUser(args.username, args.role, args.has_root)


def AddManyFastaGenomes(db, args):
    return db.AddManyFastaGenomes(
        args.batchfile, args.checkm_file, args.genome_list_id,
        args.genome_list_name, args.force
    )


def AddMarkers(db, args):
    return db.AddMarkers(
        args.batchfile, args.marker_set_id,
        args.marker_set_name, args.force
    )


def CreateTreeData(db, args):

    genome_id_list = []

    if args.all_genomes:
        genome_id_list = db.GetAllGenomeIds()
        if genome_id_list is False:
            return False
    else:
        if args.genome_ids:
            temp_list = db.ExternalGenomeIdsToGenomeIds(
                args.genome_ids.split(","))
            if temp_list is False:
                return False
            genome_id_list += temp_list

        if args.genome_list_ids:
            temp_list = db.GetGenomeIdListFromGenomeListIds(
                args.genome_list_ids.split(","))
            if temp_list is False:
                return False
            genome_id_list += temp_list

        genome_batchfile_ids = []
        if args.genome_batchfile:
            fh = open(args.genome_batchfile, "rb")
            for line in fh:
                line = line.rstrip()
                genome_batchfile_ids.append(line)

        if genome_batchfile_ids:
            genome_id_list += db.ExternalGenomeIdsToGenomeIds(
                genome_batchfile_ids)

    if (len(genome_id_list) == 0):
        db.ReportError("No genomes found from the information provided.")
        return False

    marker_id_list = []

    if args.marker_ids:
        temp_list = db.ExternalMarkerIdsToMarkerIds(args.marker_ids.split(","))
        if temp_list is None:
            return False
        marker_id_list += temp_list

    # TODO: make GetMarkerIdListFromMarkerSetIds
    if args.marker_set_ids:
        marker_set_ids = args.marker_set_ids.split(",")
        for set_id in marker_set_ids:
            temp_marker_list = db.GetMarkerIdListFromMarkerSetId(set_id)
            if temp_marker_list:
                marker_id_list += temp_marker_list

    marker_batchfile_ids = []
    if args.marker_batchfile:
        fh = open(args.marker_batchfile, "rb")
        for line in fh:
            line = line.rstrip()
            marker_batchfile_ids.append(line)

    if marker_batchfile_ids:
        temp_list = db.ExternalMarkerIdsToMarkerIds(marker_batchfile_ids)
        if temp_list is None:
            return False
        marker_id_list += temp_list

    if (len(marker_id_list) == 0):
        db.ReportError("No markers found from the information provided.")
        return False

    profile_config_dict = dict()
    if args.profile_args:
        for profile_arg in args.profile_args:
            key_value_pair = profile_arg.split('=')
            try:
                profile_config_dict[key_value_pair[0]] = key_value_pair[1]
            except IndexError:
                profile_config_dict[key_value_pair[0]] = None

    return db.MakeTreeData(marker_id_list, genome_id_list, args.out_dir, "gtdblite", args.profile, profile_config_dict, not(args.no_tree))


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

    return db.DeleteGenomes(args.batchfile, external_ids)


def CreateGenomeList(db, args):

    external_ids = []

    if args.genome_ids:
        external_ids = args.genome_ids.split(",")
    genome_list_id = db.CreateGenomeList(
        args.batchfile, external_ids, args.name, args.description, (not args.public))

    if genome_list_id is False:
        return False

    try:
        print_success = db.PrintGenomeListsDetails([genome_list_id])
    except:
        print_success = False

    if not print_success:
        db.ReportWarning(
            "New genome list was created, but failed to print details to screen.")

    return genome_list_id


def ViewGenomeLists(db, args):

    genome_lists = []

    view_all = False
    if args.view_all:
        view_all = True

    if args.root_owned or (args.self_owned and db.currentUser.isRootUser()):
        genome_lists = db.GetVisibleGenomeListsByOwner(
            all_non_private=view_all)
    elif args.self_owned:
        genome_lists = db.GetVisibleGenomeListsByOwner(
            db.currentUser.getUserId(), all_non_private=view_all)
    elif args.all_owners:
        genome_lists = db.GetAllVisibleGenomeListIds(all_non_private=view_all)
    else:
        # TODO: this
        db.ReportError(
            "Viewing other peoples' genome lists not yet implemented.")
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


def ViewMarkers(db, args):

    if args.view_all:
        return db.ViewMarkers()
    else:
        external_ids = None
        if args.id_list:
            external_ids = args.id_list.split(",")
        return db.ViewMarkers(args.batchfile, external_ids)


def CreateMarkerSet(db, args):

    external_ids = []

    if args.genome_ids:
        external_ids = args.id_list.split(",")
    marker_set_id = db.CreateMarkerSet(
        args.batchfile, external_ids, args.name, args.description, (not args.public))

    if marker_set_id is False:
        return False

    try:
        print_success = db.PrintMarkerSetsDetails([marker_set_id])
    except:
        print_success = False

    if not print_success:
        db.ReportWarning(
            "New marker set was created, but failed to print details to screen.")

    return marker_set_id


def ViewMarkerSets(db, args):

    marker_sets = []

    view_all = False
    if args.view_all:
        view_all = True

    if args.root_owned or (args.self_owned and db.currentUser.isRootUser()):
        marker_sets = db.GetVisibleMarkerSetsByOwner(all_non_private=view_all)
    elif args.self_owned:
        marker_sets = db.GetVisibleMarkerSetsByOwner(
            db.currentUser.getUserId(), all_non_private=view_all)
    elif args.all_owners:
        marker_sets = db.GetAllVisibleMarkerSetIds(all_non_private=view_all)
    else:
        # TODO: this
        db.ReportError(
            "Viewing other peoples' marker sets not yet implemented.")
        return False

    if len(marker_sets) == 0:
        print "No marker sets found."
        return True
    return db.PrintMarkerSetsDetails(marker_sets)


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


def MarkerSetsContents(db, args):
    set_ids = []

    if args.set_ids:
        set_ids = args.set_ids.split(",")

    return db.ViewMarkerSetsContents(set_ids)


def viewMetadata(db, args):
    return db.viewMetadata()


def exportMetadata(db, args):
    return db.exportMetadata(args.outfile)


def importMetadata(db, args):
    return db.importMetadata(args.table, args.field, args.typemeta, args.metadatafile)


def createMetadata(db, args):
    return db.createMetadata(args.metadatafile)

if __name__ == '__main__':

    # create the top-level parser
    parser = argparse.ArgumentParser(prog='gtdblite.py')
    parser.add_argument('-r', dest='login_as_root', action='store_true',
                        help='Login as the root user'),
    parser.add_argument('-u', dest='logon_as_user',
                        help='Logon as this user (implies -r)'),
    parser.add_argument('-t', dest='threads', type=int,
                        help='Threads to use'),
    parser.add_argument('-f', dest='force', action='store_true',
                        help='Force the action (required to override warnings for certain actions)'),
    parser.add_argument('-y', dest='assume_yes', action='store_true',
                        help='Assume yes to all confirm prompts (useful for batch processing)'),
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Run in debug mode')

    category_parser = parser.add_subparsers(
        help='Category Command Help', dest='category_parser_name')

    user_category_parser = category_parser.add_parser(
        'users', help='Access the user management sub-commands')
    user_category_subparser = user_category_parser.add_subparsers(
        help='User command help', dest='user_subparser_name')

    genome_category_parser = category_parser.add_parser(
        'genomes', help='Access the genome management sub-commands')
    genome_category_subparser = genome_category_parser.add_subparsers(
        help='Genome command help', dest='genome_subparser_name')

    genome_list_category_parser = category_parser.add_parser(
        'genome_lists', help='Access the genome list management sub-commands')
    genome_list_category_subparser = genome_list_category_parser.add_subparsers(
        help='Genome List command help', dest='genome_list_subparser_name')

    marker_category_parser = category_parser.add_parser(
        'markers', help='Access the marker management commands')
    marker_category_subparser = marker_category_parser.add_subparsers(
        help='Marker command help', dest='marker_subparser_name')

    marker_set_category_parser = category_parser.add_parser(
        'marker_sets', help='Access the marker set management sub-commands')
    marker_set_category_subparser = marker_set_category_parser.add_subparsers(
        help='Marker Set command help', dest='marker_sets_subparser_name')

    metadata_category_parser = category_parser.add_parser(
        'metadata', help='Access the metadata management commands')
    metadata_category_subparser = metadata_category_parser.add_subparsers(
        help='Metadata command help', dest='metadata_subparser_name')

    tree_category_parser = category_parser.add_parser(
        'trees', help='Access the tree management commands')
    tree_category_subparser = tree_category_parser.add_subparsers(
        help='Tree command help', dest='tree_subparser_name')

    profile_category_parser = category_parser.add_parser(
        'profiles', help='Access the profile management commands')
    profile_category_subparser = profile_category_parser.add_subparsers(
        help='Profile command help', dest='profile_subparser_name')

# -------- User Management subparsers

    # user add parser
    parser_user_add = user_category_subparser.add_parser('add',
                                                         help='Add a user')
    parser_user_add.add_argument('--username', dest='username',
                                 required=True, help='Username of the new user.')
    parser_user_add.add_argument('--role', dest='role', choices=('user', 'admin'),
                                 required=False, help='Role of the new user')
    parser_user_add.add_argument('--has_root', dest='has_root', action="store_true",
                                 required=False, help='User has permission to become the root user.')
    parser_user_add.set_defaults(func=AddUser)

    # user edit parser
    parser_user_edit = user_category_subparser.add_parser('edit',
                                                          help='Edit a user')
    parser_user_edit.add_argument('--username', dest='username',
                                  required=True, help='Username of the user to edit.')
    parser_user_edit.add_argument('--role', dest='role', choices=('user', 'admin'),
                                  required=False, help='Change the user to this role')

    mutex_group = parser_user_edit.add_mutually_exclusive_group(required=False)
    mutex_group.add_argument('--has_root', dest='has_root', action="store_true", default=None,
                             help='Grant user the permission to become the root user.')
    mutex_group.add_argument('--no_root', dest='has_root', action="store_false", default=None,
                             help="Revoke user's permission to become the root user.")
    parser_user_edit.set_defaults(func=EditUser)

    # user delete parser

# -------- Genome Management subparsers

    # genome add parser
    parser_genome_add = genome_category_subparser.add_parser('add',
                                                             help='Add one or many genomes to the tree.')
    parser_genome_add.add_argument('--batchfile', dest='batchfile',
                                   required=True, help='Batchfile describing the genomes - one genome per line, tab separated in 3-5 columns (filename, name, desc, [source], [id_at_source])')
    parser_genome_add.add_argument('--checkm_results', dest='checkm_file',
                                   required=True, help='Provide a checkM results file. MUST BE A TAB TABLE! e.g. "checkm taxonomy_wf -f CHECKM_FILE --tab_table domain Bacteria bins/ output"')

    mutex_group = parser_genome_add.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--modify_list', dest='genome_list_id',
                             help='Modify a genome list with the \
                             specified id and add all batchfile genomes into it.')
    mutex_group.add_argument('--create_list', dest='genome_list_name',
                             help='Create a genome list with the specified name and add all batchfile genomes into it.')
    mutex_group.add_argument('--no_list', dest='no_genome_list', action="store_true",
                             help="Don't add these genomes to a list.")
    parser_genome_add.set_defaults(func=AddManyFastaGenomes)

    # genome view parser
    parser_genome_view = genome_category_subparser.add_parser('view',
                                                              help='View the details of genomes in the database.')
    parser_genome_view.add_argument('--batchfile', dest='batchfile', default=None,
                                    help='Batchfile of genome ids (one per line) to view')
    parser_genome_view.add_argument('--genome_ids', dest='id_list', default=None,
                                    help='Provide a list of genome ids (comma separated) to view')
    parser_genome_view.add_argument('--all', dest='view_all', action="store_true",
                                    help='View ALL the genomes in the database. This might take a while...')
    parser_genome_view.set_defaults(func=ViewGenomes)

    # genome delete parser
    parser_genome_delete = genome_category_subparser.add_parser('delete',
                                                                help='Delete genomes in the database.')
    parser_genome_delete.add_argument('--batchfile', dest='batchfile', default=None,
                                      help='Batchfile of genome ids (one per line) to view')
    parser_genome_delete.add_argument('--genome_ids', dest='id_list', default=None,
                                      help='Provide a list of genome ids (comma separated) to view')
    parser_genome_delete.set_defaults(func=DeleteGenomes)

# -------- Genome Lists Management subparsers

    # Create genome list
    parser_genome_lists_create = genome_list_category_subparser.add_parser('create',
                                                                           help='Create a genome list')
    parser_genome_lists_create.add_argument('--batchfile', dest='batchfile',
                                            help='A file of genome IDs, one per line, to add to the create list')
    parser_genome_lists_create.add_argument('--genome_ids', dest='genome_ids',
                                            help='List of genome IDs (comma separated) to add to the create list')
    parser_genome_lists_create.add_argument('--name', dest='name', required=True,
                                            help='The name of the genome list.')
    parser_genome_lists_create.add_argument('--description', dest='description',
                                            help='A brief description of the genome list.')
    parser_genome_lists_create.add_argument('--set_public', dest='public', action='store_true', default=False,
                                            help='Make the new list publically visible.')

    parser_genome_lists_create.set_defaults(func=CreateGenomeList)

    # View genome lists
    parser_genome_lists_view = genome_list_category_subparser.add_parser('view',
                                                                         help='View visible genome lists.')

    mutex_group = parser_genome_lists_view.add_mutually_exclusive_group(
        required=True)
    mutex_group.add_argument('--root', dest='root_owned', default=False,
                             action='store_true', help='Only show genome lists owned by the root user.')
    mutex_group.add_argument('--self', dest='self_owned', default=False,
                             action='store_true', help='Only show genome lists owned by you.')
    mutex_group.add_argument(
        '--owner', dest='owner_name', help='Only show genome lists owned by a specific user.')
    mutex_group.add_argument('--everyone', dest='all_owners', default=False,
                             action='store_true', help='Show genome lists of all users.')

    parser_genome_lists_view.add_argument('--all_accessible', dest='view_all', default=False, action='store_true',
                                          help='View EVERY genome list that you can access (the default is to show only public and your own lists).')

    parser_genome_lists_view.set_defaults(func=ViewGenomeLists)

    # Show genome list
    parser_genome_lists_contents = genome_list_category_subparser.add_parser('contents',
                                                                             help='View the contents of genome list(s)')
    parser_genome_lists_contents.add_argument('--list_ids', dest='list_ids', required=True,
                                              help='Provide a list of genome list ids (comma separated) whose contents you wish to view.')
    parser_genome_lists_contents.set_defaults(func=ContentsGenomeLists)

    # Edit genome list
    parser_genome_lists_edit = genome_list_category_subparser.add_parser('edit',
                                                                         help='Edit a genome list')
    parser_genome_lists_edit.add_argument('--list_id', dest='list_id',
                                          required=True, help='The id of the genome list to edit')
    parser_genome_lists_edit.add_argument('--batchfile', dest='batchfile',
                                          help='A file of genome IDs, one per line, to add remove from the list')
    parser_genome_lists_edit.add_argument('--genome_ids', dest='genome_ids',
                                          help='List of genome IDs to add/remove from list')
    parser_genome_lists_edit.add_argument('--operation', dest='operation', choices=('add', 'remove'),
                                          help='What to do with the tree_ids with regards to the genome list.')
    parser_genome_lists_edit.add_argument('--name', dest='name',
                                          help='Modify the name of the list to this.')
    parser_genome_lists_edit.add_argument('--description', dest='description',
                                          help='Change the brief description of the genome list to this.')

    mutex_group = parser_genome_lists_edit.add_mutually_exclusive_group(
        required=False)
    mutex_group.add_argument('--set_private', dest='private', action="store_true", default=False,
                             help='Make this genome list private (only you can see).')
    mutex_group.add_argument('--set_public', dest='public', action="store_true", default=False,
                             help='Make this genome list public (all users can see).')
    parser_genome_lists_edit.set_defaults(func=EditGenomeLists)

# --------- Marker Management Subparsers

    parser_marker_add = marker_category_subparser.add_parser('add',
                                                             help='Add in one or many marker HMMs into the database')
    parser_marker_add.add_argument('--batchfile', dest='batchfile', required=True,
                                   help='Batchfile describing the markers - one HMM file per line (one model per file), tab separated in 3-5 columns (filename, name, desc, [database], [database_specific_id]')
    mutex_group = parser_marker_add.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--modify_set', dest='marker_set_id',
                             help='Modify a marker set with the \
                            specified id and add all markers to it.')
    mutex_group.add_argument('--create_set', dest='marker_set_name',
                             help='Create a marker set with the specified name and add these markers to it.')
    mutex_group.add_argument('--no_set', dest='no_marker_set', action="store_true",
                             help="Don't add these markers to a marker set.")
    parser_marker_add.set_defaults(func=AddMarkers)

    # View markers
    parser_marker_view = marker_category_subparser.add_parser('view',
                                                              help='View HMM markers in the database')
    parser_marker_view.add_argument('--batchfile', dest='batchfile', default=None,
                                    help='Batchfile of marker ids (one per line) to view')
    parser_marker_view.add_argument('--marker_ids', dest='id_list', default=None,
                                    help='Provide a list of genome ids (comma separated) to view')
    parser_marker_view.add_argument('--all', dest='view_all', action="store_true",
                                    help='View ALL the markers in the database.')
    parser_marker_view.set_defaults(func=ViewMarkers)

# -------- Marker Set Management subparsers

    # Create marker set
    parser_marker_sets_create = marker_set_category_subparser.add_parser('create',
                                                                         help='Create a marker set')
    parser_marker_sets_create.add_argument('--batchfile', dest='batchfile',
                                           help='A file of marker IDs, one per line, to add to the created set')
    parser_marker_sets_create.add_argument('--marker_ids', dest='genome_ids',
                                           help='List of marker IDs (comma separated) to add to the created set')
    parser_marker_sets_create.add_argument('--name', dest='name', required=True,
                                           help='The name of the marker set.')
    parser_marker_sets_create.add_argument('--description', dest='description',
                                           help='A brief description of the marker set.')
    parser_marker_sets_create.add_argument('--set_public', dest='public', action='store_true', default=False,
                                           help='Make the new set publically visible.')

    parser_marker_sets_create.set_defaults(func=CreateMarkerSet)

    parser_marker_sets_view = marker_set_category_subparser.add_parser('view',
                                                                       help='View visible marker sets.')

    mutex_group = parser_marker_sets_view.add_mutually_exclusive_group(
        required=True)
    mutex_group.add_argument('--root', dest='root_owned', default=False,
                             action='store_true', help='Only show marker sets owned by the root user.')
    mutex_group.add_argument('--self', dest='self_owned', default=False,
                             action='store_true', help='Only show marker sets owned by you.')
    mutex_group.add_argument(
        '--owner', dest='owner_name', help='Only show marker sets owned by a specific user.')
    mutex_group.add_argument('--everyone', dest='all_owners', default=False,
                             action='store_true', help='Show marker sets of all users.')

    parser_marker_sets_view.add_argument('--all_accessible', dest='view_all', default=False, action='store_true',
                                         help='View EVERY marker set that you can access (the default is to show only public and your own sets).')

    parser_marker_sets_view.set_defaults(func=ViewMarkerSets)

    # Show marker set(s) contents
    parser_marker_sets_contents = marker_set_category_subparser.add_parser('contents',
                                                                           help='View the contents of marker set(s)')
    parser_marker_sets_contents.add_argument('--set_ids', dest='set_ids', required=True,
                                             help='Provide a list of marker set ids (comma separated) whose contents you wish to view.')
    parser_marker_sets_contents.set_defaults(func=MarkerSetsContents)

    # Edit marker set
    parser_marker_sets_edit = marker_set_category_subparser.add_parser('edit',
                                                                       help='Edit a marker set')
    parser_marker_sets_edit.add_argument('--set_id', dest='set_id',
                                         required=True, help='The id of the marker set to edit')
    parser_marker_sets_edit.add_argument('--batchfile', dest='batchfile',
                                         help='A file of marker ids, one per line, to add remove from the set')
    parser_marker_sets_edit.add_argument('--marker_ids', dest='marker_ids',
                                         help='List (comma separated) of marker ids to add/remove from set')
    parser_marker_sets_edit.add_argument('--operation', dest='operation', choices=('add', 'remove'),
                                         help='What to do with the provided marker ids with regards to the marker set.')
    parser_marker_sets_edit.add_argument('--name', dest='name',
                                         help='Modify the name of the set to this.')
    parser_marker_sets_edit.add_argument('--description', dest='description',
                                         help='Change the brief description of the marker set to this.')

    mutex_group = parser_marker_sets_edit.add_mutually_exclusive_group(
        required=False)
    mutex_group.add_argument('--set_private', dest='private', action="store_true", default=False,
                             help='Make this marker set private (only you can see).')
    mutex_group.add_argument('--set_public', dest='public', action="store_true", default=False,
                             help='Make this marker set public (all users can see).')

    parser_marker_sets_edit.set_defaults(func=EditMarkerSet)

# -------- Metadata Management subparsers

    # metadata create columns parser
    parser_metadata_create = metadata_category_subparser.add_parser('create',
                                                                    help='Create one or many metadata field.')
    parser_metadata_create.add_argument('--file', dest='metadatafile',
                                        required=True, help='Metadata file describing the new field \
                                        - one metadata per line, tab separated in 4 columns (name,description,datatype,metadata table')
    parser_metadata_create.set_defaults(func=createMetadata)

    # metadata view parser
    parser_metadata_view = metadata_category_subparser.add_parser('view',
                                                                  help='List all existing metadata fields with table name and description')
    parser_metadata_view.set_defaults(func=viewMetadata)

    # metadata import parser
    parser_metadata_import = metadata_category_subparser.add_parser('import',
                                                                    help='Import Metadata values for a list of genome')
    parser_metadata_import.add_argument('--table', dest='table', default=None, required=True,
                                        help='Table where the metadata field is present')
    parser_metadata_import.add_argument('--field', dest='field', default=None, required=True,
                                        help='Metadata field where the value(s) will be saved')
    parser_metadata_import.add_argument('--type', dest='typemeta', default=None, required=True,
                                        help='Type of the Metadata field')
    parser_metadata_import.add_argument('--metadatafile', dest='metadatafile', default=None,
                                        help='TSV file . One genome per line , tab separated in 2 columns (genome id , metadata value)')
    parser_metadata_import.set_defaults(func=importMetadata)

    # metadata export parser

    parser_metadata_export = metadata_category_subparser.add_parser('export',
                                                                    help='Export a TSV file with all Metadata fields')
    parser_metadata_export.add_argument('--output', dest='outfile', default=None, required=True,
                                        help='Destination to write the TSV file')
    parser_metadata_export.set_defaults(func=exportMetadata)

# -------- Generate Tree Data
    parser_tree_create = tree_category_subparser.add_parser('create',
                                                            help='Create a genome tree')

    parser_tree_create.add_argument('--genome_batchfile', dest='genome_batchfile', default=None,
                                    help='Provide a file of genome IDs, one per line, to include in the tree')
    parser_tree_create.add_argument('--genome_ids', dest='genome_ids', default=None,
                                    help='Provide a list of genome ids (comma separated), whose genomes should be included in the tree.')
    parser_tree_create.add_argument('--genome_list_ids', dest='genome_list_ids', default=None,
                                    help='Provide a list of genome list ids (comma separated), whose genomes should be included in the tree.')
    parser_tree_create.add_argument('--all_genomes', dest='all_genomes', default=False, action='store_true',
                                    help='Included ALL genomes in the database in the created tree.')

    parser_tree_create.add_argument('--marker_batchfile', dest='marker_batchfile', default=None,
                                    help='Provide a file of marker IDs, one per line, to build the tree')
    parser_tree_create.add_argument('--marker_ids', dest='marker_ids', default=None,
                                    help='Provide a list of marker ids (comma separated), whose markers will be used to build the tree.')
    parser_tree_create.add_argument('--marker_set_ids', dest='marker_set_ids', default=None,
                                    help='Provide a list of marker set ids (comma separated), whose markers will be used to build the tree.')

    parser_tree_create.add_argument('--output', dest='out_dir', required=True,
                                    help='Directory to output the files')
    parser_tree_create.add_argument('--no_tree', dest='no_tree', action="store_true",
                                    help="Only output tree data, don't build the tree.")

    parser_tree_create.add_argument('--profile', dest='profile',
                                    help='Tree creation profile to use (default: %s)' % (profiles.ReturnDefaultProfileName(),))
    parser_tree_create.add_argument('--profile_args', dest='profile_args', nargs='+',
                                    help='Arguments to provide to the profile')

    parser_tree_create.set_defaults(func=CreateTreeData)

    # Do the parsing
    args = parser.parse_args()

    # Special parser checks
    if (args.category_parser_name == 'trees' and args.tree_subparser_name == 'create'):
        if (args.genome_batchfile is None and args.genome_ids is None and args.genome_list_ids is None and not args.all_genomes):
            parser_tree_create.error(
                'Need to specify at least one of --genome_batchfile, --genome_ids, --genome_list_ids or --all_genomes')
        if (args.marker_batchfile is None and args.marker_ids is None and args.marker_set_ids is None):
            parser_tree_create.error(
                'Need to specify at least one of --marker_batchfile, --marker_ids or --marker_set_ids')

    if (args.category_parser_name == 'genomes' and args.genome_subparser_name == 'view'):
        if (args.batchfile is not None or args.id_list is not None):
            if args.view_all:
                parser_genome_view.error(
                    'argument --all must be used by itself')
        elif not args.view_all:
            parser_genome_view.error(
                'need to specify at least one of --all, --batchfile or --genome_ids')

    if (args.category_parser_name == 'genomes' and args.genome_subparser_name == 'delete'):
        if (args.batchfile is None and args.id_list is None):
            parser_genome_delete.error(
                'need to specify at least one of --batchfile or --genome_ids')

    if (args.category_parser_name == 'markers' and args.marker_subparser_name == 'view'):
        if (args.batchfile is not None or args.id_list is not None):
            if args.view_all:
                parser_marker_view.error(
                    'argument --all must be used by itself')
        elif not args.view_all:
            parser_marker_view.error(
                'need to specify at least one of --all, --batchfile or --marker_ids')

    # Initialise the backend
    if args.threads:
        db = GenomeDatabase.GenomeDatabase(args.threads)
    else:
        db = GenomeDatabase.GenomeDatabase()

    db.conn.MakePostgresConnection()

    if args.debug:
        db.SetDebugMode(True)

    # Login
    try:
        if args.logon_as_user:
            if not db.RootLogin(GetLinuxUsername()):
                raise GenomeDatabaseError(
                    "Unable to impersonate user %s." % args.logon_as_user)

            if not db.UserLogin(args.logon_as_user):
                raise GenomeDatabaseError(
                    "Unable to impersonate user %s." % args.logon_as_user)

        elif args.login_as_root:
            if not db.RootLogin(GetLinuxUsername()):
                raise GenomeDatabaseError("Unable to become root user.")
        else:
            if not db.UserLogin(GetLinuxUsername()):
                raise GenomeDatabaseError("Database login failed.")

    except GenomeDatabaseError as e:
        ErrorReport(e.message + " The following error(s) were reported:\n")
        DumpDBErrors(db)
        sys.exit(-1)

    try:
        result = args.func(db, args)
    except:
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
        sys.exit(-1)
