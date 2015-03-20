#!/usr/bin/env python
import argparse
import sys
import os
import pwd

from gtdblite import GenomeDatabase
from gtdblite import profiles

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

    if args.genome_ids:
        temp_list = db.ExternalGenomeIdsToGenomeIds(args.genome_ids.split(","))
        if temp_list is None:
            return False
        genome_id_list += temp_list

    if args.genome_list_ids:
        temp_list = db.GetGenomeIdListFromGenomeListIds(args.genome_list_ids.split(","))
        if temp_list is None:
            return False
        genome_id_list += temp_list
    
    batchfile_ids = []
    if args.genome_batchfile:
        fh = open(args.batchfile, "rb")
        for line in fh:
            line = line.rstrip()
            batchfile_ids.append(line)

    if batchfile_ids:
        temp_list = db.GetGenomeIdListFromGenomeListIds(args.genome_list_ids.split(","))
        genome_id_list += db.ExternalGenomeIdsToGenomeIds(batchfile_ids)

    if (len(genome_id_list) == 0):
        db.ReportError("No genomes found from the information provided.")
        return False

    marker_id_list = []

    if args.marker_ids:
        marker_id_list += db.ExternalMarkerIdsToMarkerIds(args.marker_ids.split(","))

    # TODO: make GetMarkerIdListFromMarkerSetIds
    if args.marker_set_ids:
        marker_set_ids = args.marker_set_ids.split(",")
        for set_id in marker_set_ids:
            temp_marker_list = db.GetMarkerIdListFromMarkerSetId(set_id)
            if temp_marker_list:
                marker_id_list += temp_marker_list

    batchfile_ids = []
    if args.marker_batchfile:
        fh = open(args.batchfile, "rb")
        for line in fh:
            line = line.rstrip()
            batchfile_ids.append(line)

    if batchfile_ids:
        marker_id_list += db.ExternalMarkerIdsToMarkerIds(batchfile_ids)

    if (len(marker_id_list) == 0):
        db.ReportError("No markers found from the information provided.")
        return False

    profile_config_dict = dict()
    if args.profile_args:
        profile_args = args.profile_args.split(',')
        for profile_arg in profile_args:
            key_value_pair = profile_arg.split('=')
            try:
                profile_config_dict[key_value_pair[0]] = key_value_pair[1]
            except IndexError:
                profile_config_dict[key_value_pair[0]] = None

    
    return db.MakeTreeData(marker_id_list, genome_id_list, args.out_dir, "prefix", args.profile, profile_config_dict, not(args.no_tree))

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


def ViewGenomeLists(db, args):

    genome_lists = []

    if args.root_owned or (args.self_owned and db.currentUser.isRootUser()):
        genome_lists = db.GetVisibleGenomeListsByOwner()
    elif args.self_owned:
        genome_lists = db.GetVisibleGenomeListsByOwner(db.currentUser.getUserId())
    elif args.show_all:
        genome_lists = db.GetAllVisibleGenomeListIds()
    else:
        # TODO: this
        db.ReportError("Viewing other peoples' genome lists not yet implemented.")
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

    
def ViewMarkerSets(db, args):

    marker_sets = []

    if args.root_owned or (args.self_owned and db.currentUser.isRootUser()):
        marker_sets = db.GetVisibleMarkerSetsByOwner()
    elif args.self_owned:
        marker_sets = db.GetVisibleMarkerSetsByOwner(db.currentUser.getUserId())
    elif args.show_all:
        marker_sets = db.GetAllVisibleMarkerSetIds()
    else:
        # TODO: this
        db.ReportError("Viewing other peoples' marker sets not yet implemented.")
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

if __name__ == '__main__':

    # create the top-level parser
    parser = argparse.ArgumentParser(prog='gtdblite.py')
    parser.add_argument('-r', dest='login_as_root', action='store_true',
                        help='Login as the root user'),
    parser.add_argument('-f', dest='force', action='store_true',
                        help='Force the action (required to override warnings for certain actions)'),
    parser.add_argument('-y', dest='assume_yes', action='store_true',
                        help='Assume yes to all confirm prompts (useful for batch processing)'),
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Run in debug mode')

    category_parser = parser.add_subparsers(help='Category Command Help', dest='category_parser_name')

    user_category_parser = category_parser.add_parser('users', help='Access the user management sub-commands')
    user_category_subparser = user_category_parser.add_subparsers(help='User command help', dest='user_subparser_name')

    genome_category_parser = category_parser.add_parser('genomes', help='Access the genome management sub-commands')
    genome_category_subparser = genome_category_parser.add_subparsers(help='Genome command help', dest='genome_subparser_name')

    genome_list_category_parser = category_parser.add_parser('genome_lists', help='Access the genome list management sub-commands')
    genome_list_category_subparser = genome_list_category_parser.add_subparsers(help='Genome List command help', dest='genome_list_subparser_name')

    marker_category_parser = category_parser.add_parser('markers', help='Access the marker management commands')
    marker_category_subparser = marker_category_parser.add_subparsers(help='Marker command help', dest='marker_subparser_name')

    marker_set_category_parser = category_parser.add_parser('marker_sets', help='Access the marker set management sub-commands')
    marker_set_category_subparser = marker_set_category_parser.add_subparsers(help='Marker Set command help', dest='marker_sets_subparser_name')

    tree_category_parser = category_parser.add_parser('trees', help='Access the tree management commands')
    tree_category_subparser = tree_category_parser.add_subparsers(help='Tree command help', dest='tree_subparser_name')

    profile_category_parser = category_parser.add_parser('profiles', help='Access the profile management commands')
    profile_category_subparser = profile_category_parser.add_subparsers(help='Profile command help', dest='profile_subparser_name')

# -------- User Management subparsers

    # user add parser
    parser_user_add = user_category_subparser.add_parser('add',
                                    help='Add a user')
    parser_user_add.add_argument('--username', dest = 'username',
                                    required=True, help='Username of the new user.')
    parser_user_add.add_argument('--role', dest = 'role', choices = ('user', 'admin'), 
                                    required=False, help='Role of the new user')
    parser_user_add.add_argument('--has_root', dest = 'has_root', action="store_true", 
                                    required=False, help='User has permission to become the root user.')
    parser_user_add.set_defaults(func=AddUser)
    
    
    # user edit parser
    parser_user_edit = user_category_subparser.add_parser('edit',
                                    help='Edit a user')
    parser_user_edit.add_argument('--username', dest = 'username',
                                    required=True, help='Username of the user to edit.')
    parser_user_edit.add_argument('--role', dest = 'role', choices = ('user', 'admin'), 
                                    required=False, help='Change the user to this role')
    
    mutex_group = parser_user_edit.add_mutually_exclusive_group(required=False)
    mutex_group.add_argument('--has_root', dest = 'has_root', action="store_true", default=None,
                                    help='Grant user the permission to become the root user.')
    mutex_group.add_argument('--no_root', dest = 'has_root', action="store_false", default=None,
                                    help="Revoke user's permission to become the root user.")
    parser_user_edit.set_defaults(func=EditUser)
    
    
    # user delete parser
    
# -------- Genome Management subparsers

    # genome add parser
    parser_genome_add = genome_category_subparser.add_parser('add',
                                    help='Add one or many genomes to the tree.')
    parser_genome_add.add_argument('--batchfile', dest = 'batchfile',
                                    required=True, help='Batchfile describing the genomes - one genome per line, tab separated in 3-5 columns (filename, name, desc, [source], [id_at_source])')
    parser_genome_add.add_argument('--checkm_results', dest = 'checkm_file',
                                    required=True, help='Provide a checkM results file. MUST BE A TAB TABLE! e.g. "checkm taxonomy_wf -f CHECKM_FILE --tab_table domain Bacteria bins/ output"')

    mutex_group = parser_genome_add.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--modify_list', dest = 'genome_list_id',
                                    help='Modify a genome list with the \
                                    specified id and add all batchfile genomes into it.')
    mutex_group.add_argument('--create_list', dest = 'genome_list_name',
                                    help='Create a genome list with the specified name and add all batchfile genomes into it.')
    mutex_group.add_argument('--no_list', dest = 'no_genome_list', action="store_true",
                                    help="Don't add these genomes to a list.")
    parser_genome_add.set_defaults(func=AddManyFastaGenomes)

    # genome view parser
    parser_genome_view = genome_category_subparser.add_parser('view',
                                    help='View the details of genomes in the database.')
    parser_genome_view.add_argument('--batchfile', dest = 'batchfile', default=None,
                                    help='Batchfile of genome ids (one per line) to view')
    parser_genome_view.add_argument('--genome_ids', dest = 'id_list', default=None,
                                    help='Provide a list of genome ids (comma separated) to view')
    parser_genome_view.add_argument('--all', dest = 'view_all', action="store_true",
                                    help='View ALL the genomes in the database. This might take a while...')
    parser_genome_view.set_defaults(func=ViewGenomes)

    # genome delete parser
    parser_genome_delete = genome_category_subparser.add_parser('delete',
                                    help='Delete genomes in the database.')
    parser_genome_delete.add_argument('--batchfile', dest = 'batchfile', default=None,
                                    help='Batchfile of genome ids (one per line) to view')
    parser_genome_delete.add_argument('--genome_ids', dest = 'id_list', default=None,
                                    help='Provide a list of genome ids (comma separated) to view')
    parser_genome_delete.set_defaults(func=DeleteGenomes)

# -------- Genome Lists Management subparsers

    #------------ View genome lists
    parser_genome_lists_view = genome_list_category_subparser.add_parser('view',
                                        help='View visible genome lists.')

    mutex_group = parser_genome_lists_view.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--root', dest = 'root_owned', default=False,
                                        action='store_true', help='Only show genome lists owned by the root user.')
    mutex_group.add_argument('--self', dest = 'self_owned', default=False,
                                        action='store_true', help='Only show genome lists owned by you.')
    mutex_group.add_argument('--owner', dest = 'owner_name', help='Only show genome lists owned by a specific user.')
    mutex_group.add_argument('--all', dest = 'show_all', default=False,
                             action='store_true', help='Show all visible genome lists')

    parser_genome_lists_view.set_defaults(func=ViewGenomeLists)

    #------------ Show genome list
    parser_genome_lists_contents = genome_list_category_subparser.add_parser('contents',
                                        help='View the contents of genome list(s)')
    parser_genome_lists_contents.add_argument('--list_ids', dest = 'list_ids', required=True,
                                        help='Provide a list of genome list ids (comma separated) whose contents you wish to view.')
    parser_genome_lists_contents.set_defaults(func=ContentsGenomeLists)

    #------------ Edit genome list
    parser_genome_lists_edit = genome_list_category_subparser.add_parser('edit',
                                        help='Edit a genome list') 
    parser_genome_lists_edit.add_argument('--list_id', dest = 'list_id',
                                        required=True, help='The id of the genome list to edit')
    parser_genome_lists_edit.add_argument('--batchfile', dest = 'batchfile',
                                        help='A file of genome IDs, one per line, to add remove from the list')
    parser_genome_lists_edit.add_argument('--genome_ids', dest = 'genome_ids',
                                        help='List of tree_ids to add/remove from list')
    parser_genome_lists_edit.add_argument('--operation', dest = 'operation', choices=('add','remove'),
                                        help='What to do with the tree_ids with regards to the genome list.')
    parser_genome_lists_edit.add_argument('--name', dest = 'name',
                                        help='Modify the name of the list to this.')
    parser_genome_lists_edit.add_argument('--description', dest = 'description',
                                        help='Change the brief description of the genome list to this.')
    
    mutex_group = parser_genome_lists_edit.add_mutually_exclusive_group(required=False)
    mutex_group.add_argument('--set_private', dest = 'private', action="store_true", default=False,
                             help='Make this genome list private (only you can see).')
    mutex_group.add_argument('--set_public', dest = 'public', action="store_true", default=False,
                             help='Make this genome list public (all users can see).')
    parser_genome_lists_edit.set_defaults(func=EditGenomeLists)

#--------- Marker Management Subparsers

    parser_marker_add = marker_category_subparser.add_parser('add',
                                 help='Add in one or many marker HMMs into the database')
    parser_marker_add.add_argument('--batchfile', dest='batchfile', required=True,
                                help='Batchfile describing the markers - one HMM file per line (one model per file), tab separated in 3-5 columns (filename, name, desc, [database], [database_specific_id]')
    mutex_group = parser_marker_add.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--modify_set', dest = 'marker_set_id',
                                    help='Modify a marker set with the \
                                    specified id and add all markers to it.')
    mutex_group.add_argument('--create_set', dest = 'marker_set_name',
                                    help='Create a marker set with the specified name and add these markers to it.')
    mutex_group.add_argument('--no_set', dest = 'no_marker_set', action="store_true",
                                    help="Don't add these markers to a marker set.")
    parser_marker_add.set_defaults(func=AddMarkers)

    # View markers
    parser_marker_view = marker_category_subparser.add_parser('view',
                                 help='View HMM markers in the database')
    parser_marker_view.add_argument('--batchfile', dest = 'batchfile', default=None,
                                    help='Batchfile of marker ids (one per line) to view')
    parser_marker_view.add_argument('--marker_ids', dest = 'id_list', default=None,
                                    help='Provide a list of genome ids (comma separated) to view')
    parser_marker_view.add_argument('--all', dest = 'view_all', action="store_true",
                                    help='View ALL the markers in the database.')
    parser_marker_view.set_defaults(func=ViewMarkers)

# -------- Marker Set Management subparsers

    parser_marker_sets_view = marker_set_category_subparser.add_parser('view',
                                        help='View visible marker sets.')

    mutex_group = parser_marker_sets_view.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--root', dest = 'root_owned', default=False,
                                        action='store_true', help='Only show marker sets owned by the root user.')
    mutex_group.add_argument('--self', dest = 'self_owned', default=False,
                                        action='store_true', help='Only show marker sets owned by you.')
    mutex_group.add_argument('--owner', dest = 'owner_name', help='Only show marker sets owned by a specific user.')
    mutex_group.add_argument('--all', dest = 'show_all', default=False,
                             action='store_true', help='Show all marker sets.')

    parser_marker_sets_view.set_defaults(func=ViewMarkerSets)


    #------------ Show marker set(s) contents
    parser_marker_sets_contents = marker_set_category_subparser.add_parser('contents',
                                        help='View the contents of marker set(s)')
    parser_marker_sets_contents.add_argument('--set_ids', dest = 'set_ids', required=True,
                                        help='Provide a list of marker set ids (comma separated) whose contents you wish to view.')
    parser_marker_sets_contents.set_defaults(func=MarkerSetsContents)

    #------------ Edit marker set 
    parser_marker_sets_edit = marker_set_category_subparser.add_parser('edit',
                                        help='Edit a marker set') 
    parser_marker_sets_edit.add_argument('--set_id', dest = 'set_id',
                                        required=True, help='The id of the marker set to edit')
    parser_marker_sets_edit.add_argument('--batchfile', dest = 'batchfile',
                                        help='A file of marker ids, one per line, to add remove from the set')
    parser_marker_sets_edit.add_argument('--marker_ids', dest = 'marker_ids',
                                        help='List (comma separated) of marker ids to add/remove from set')
    parser_marker_sets_edit.add_argument('--operation', dest = 'operation', choices=('add','remove'),
                                        help='What to do with the provided marker ids with regards to the marker set.')
    parser_marker_sets_edit.add_argument('--name', dest = 'name',
                                        help='Modify the name of the set to this.')
    parser_marker_sets_edit.add_argument('--description', dest = 'description',
                                        help='Change the brief description of the marker set to this.')
    
    mutex_group = parser_marker_sets_edit.add_mutually_exclusive_group(required=False)
    mutex_group.add_argument('--set_private', dest = 'private', action="store_true", default=False,
                             help='Make this marker set private (only you can see).')
    mutex_group.add_argument('--set_public', dest = 'public', action="store_true", default=False,
                             help='Make this marker set public (all users can see).')
    
    parser_marker_sets_edit.set_defaults(func=EditMarkerSet)

# -------- Generate Tree Data

    parser_tree_create = tree_category_subparser.add_parser('create',
                                        help='Create a genome tree')

    parser_tree_create.add_argument('--genome_batchfile', dest = 'genome_batchfile', default=None,
                                        help='Provide a file of genome IDs, one per line, to include in the tree')
    parser_tree_create.add_argument('--genome_ids', dest = 'genome_ids', default=None,
                                        help='Provide a list of genome ids (comma separated), whose genomes should be included in the tree.')
    parser_tree_create.add_argument('--genome_list_ids', dest = 'genome_list_ids', default=None,
                                        help='Provide a list of genome list ids (comma separated), whose genomes should be included in the tree.')

    parser_tree_create.add_argument('--marker_batchfile', dest = 'marker_batchfile', default=None,
                                        help='Provide a file of marker IDs, one per line, to build the tree')
    parser_tree_create.add_argument('--marker_ids', dest = 'marker_ids', default=None,
                                        help='Provide a list of marker ids (comma separated), whose markers will be used to build the tree.')
    parser_tree_create.add_argument('--marker_set_ids', dest = 'marker_set_ids', default=None,
                                        help='Provide a list of marker set ids (comma separated), whose markers will be used to build the tree.')

    parser_tree_create.add_argument('--output', dest = 'out_dir', required=True,
                                        help='Directory to output the files')
    parser_tree_create.add_argument('--no_tree', dest = 'no_tree', action="store_true",
                                        help="Only output tree data, don't build the tree.")

    parser_tree_create.add_argument('--profile', dest = 'profile',
                                        help='Tree creation profile to use (default: %s)' % (profiles.ReturnDefaultProfileName(),))
    parser_tree_create.add_argument('--profile_args', dest = 'profile_args',
                                        help='Arguments to provide to the profile')

    parser_tree_create.set_defaults(func=CreateTreeData)

    """
# -------- Create users

    parser_createuser = subparsers.add_parser('CreateUser',
                                              help='Create user help')
    parser_createuser.add_argument('--user', dest = 'username',
                                   required=True, help='Username of the created user')
    parser_createuser.add_argument('--type', dest = 'type',
                                   required=True, help='User type')
    parser_createuser.set_defaults(func=CreateUser)

# -------- Modify users

    parser_modifyuser = subparsers.add_parser('ModifyUser', help='Modify user help')
    parser_modifyuser.add_argument('--user', dest = 'username',
                                   required=True, help='Username of the user')
    parser_modifyuser.add_argument('--type', dest = 'type', help='User type')
    parser_modifyuser.add_argument('--password', dest = 'password',
                                   action = 'store_true', help='User type')
    parser_modifyuser.set_defaults(func=ModifyUser)

# -------- Show users

    parser_showuser = subparsers.add_parser('ShowUser', help='Show user help')
    parser_showuser.add_argument('--user', dest = 'username',
                                required=True, help='Username of the user')
    parser_showuser.set_defaults(func=ShowUser)

# -------- Delete users

    parser_deleteuser = subparsers.add_parser('DeleteUser', help='Delete user help')
    parser_deleteuser.add_argument('--user', dest = 'username',
                                   required=True, help='Username of the user to delete')
    parser_deleteuser.add_argument('--force', dest = 'force', action='store_true',
                                   help='Do not prompt for confirmation')
    parser_deleteuser.set_defaults(func=DeleteUser)

# -------- Genome management subparsers

    parser_addfastagenome = subparsers.add_parser('AddFastaGenome',
                                    help='Add a genome to the tree from a Fasta file')
    parser_addfastagenome.add_argument('--file', dest = 'filename',
                                       required=True, help='FASTA file to add')
    parser_addfastagenome.add_argument('--name', dest = 'name',
                                       required=True, help='Name of the genome')
    parser_addfastagenome.add_argument('--description', dest = 'description',
                                       required=True, help='Brief description of the genome')
    parser_addfastagenome.add_argument('--source', dest = 'source',
                                       help='The source of this genome (see ShowGenomeSources)')
    parser_addfastagenome.add_argument('--modify_list', dest = 'genome_list_id',
                                       help='Modify a genome list with the \
                                       specified id by adding the current \
                                       genome to it')
    parser_addfastagenome.add_argument('--id_at_source', dest = 'id_at_source',
                                       help='The id of this genome at the specified source')
    parser_addfastagenome.set_defaults(func=AddFastaGenome)


    parser_addmanyfastagenomes = subparsers.add_parser('AddManyFastaGenomes',
                                    help='Add a genome to the tree from a Fasta file')
    parser_addmanyfastagenomes.add_argument('--batchfile', dest = 'batchfile',
                                    required=True, help='Add genomes en masse with a batch file (one genome per line, tab separated in 3 columns (filename,name,desc))')
    mutex_group = parser_addmanyfastagenomes.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--modify_list', dest = 'genome_list_id',
                                    help='Modify a genome list with the \
                                    specified id and add all batchfile genomes into it.')
    mutex_group.add_argument('--create_list', dest = 'genome_list_name',
                                    help='Create a genome list with the specified name and add all batchfile genomes into it.')
    parser_addmanyfastagenomes.set_defaults(func=AddManyFastaGenomes)

# --------- Export FASTA Genome

    parser_exportfasta = subparsers.add_parser('ExportFasta',
                                    help='Export a genome to a FASTA file')
    parser_exportfasta.add_argument('--tree_id', dest = 'tree_id', action='append',
                                    help='Tree ID to export. This option can be '
                                    'specified multiple times')
    parser_exportfasta.add_argument('--outdir', dest = 'prefix', default='.',
                                    help='output directory to use when exporting genomes with a batch file')
    parser_exportfasta.add_argument('--output', dest = 'output_fasta',
                                    help='Output the genome to a FASTA file')
    parser_exportfasta.add_argument('--batchfile', dest='batchfile',
                                    help='A file containing tree ids to extract')
    parser_exportfasta.set_defaults(func=ExportFasta)


# --------- Delete FASTA Genome

    parser_deletegenome = subparsers.add_parser('DeleteGenome',
                                    help='Delete a genome from the database')
    parser_deletegenome.add_argument('--tree_ids', dest = 'tree_ids',
                                    required=True, help='List of Tree IDs (comma separated)')
    parser_deletegenome.set_defaults(func=DeleteGenome)

# --------- Genome Searching

    parser_searchgenome = subparsers.add_parser('SearchGenomes',
                                    help='Add a genome to the tree from a Fasta file')
    parser_searchgenome.add_argument('--name', dest = 'name',
                                       help='Search for genomes containing this name')
    parser_searchgenome.add_argument('--description', dest = 'description',
                                       help='Search for genomes containing this description')
    parser_searchgenome.add_argument('--tree_id', dest = 'tree_id',
                                       help='Show genome with this tree_id')
    parser_searchgenome.add_argument('--list_id', dest = 'list_id',
                                       help='Show all genomes in this list')
    parser_searchgenome.add_argument('--owner', dest = 'owner', nargs='?', default='-1',
                                       help='Search for genomes owned by this username. ' +
                                      'With no parameter finds genomes owned by the current user')
    parser_searchgenome.set_defaults(func=SearchGenomes)

# --------- Show Genome Sources

    parser_showgenomesources = subparsers.add_parser('ShowGenomeSources',
                                help='Show the sources of the genomes')
    parser_showgenomesources.set_defaults(func=ShowGenomeSources)

# --------- Create A Genome List

    parser_creategenomelist = subparsers.add_parser('CreateGenomeList',
                                        help='Create a genome list from a list of accessions')
    parser_creategenomelist.add_argument('--file', dest = 'filename',
                                       required=True, help='File containing list of accessions')
    parser_creategenomelist.add_argument('--source', dest = 'source',
                                       help='Source of the accessions listed in the file')
    parser_creategenomelist.add_argument('--name', dest = 'name',
                                       required=True, help='Name of the genome list')
    parser_creategenomelist.add_argument('--description', dest = 'description',
                                       required=True, help='Brief description of the genome list')
    parser_creategenomelist.add_argument('--public', dest = 'public', default=False,
                                       action='store_true', help='Make the list visible to all users.')
    parser_creategenomelist.set_defaults(func=CreateGenomeList)

# --------- Modify A Genome List

    parser_modifygenomelist = subparsers.add_parser('ModifyGenomeList',
                                        help='Modify a genome list')
    parser_modifygenomelist.add_argument('--list_id', dest = 'list_id',
                                        required=True, help='File containing list of accessions')
    parser_modifygenomelist.add_argument('--tree_ids', dest = 'tree_ids',
                                        help='List of tree_ids to add/remove from list')
    parser_modifygenomelist.add_argument('--operation', dest = 'operation', choices=('add','remove'),
                                        help='What to do with the tree_ids with regards to the genome list.')
    parser_modifygenomelist.add_argument('--description', dest = 'description',
                                        help='Change the brief description of the genome list to this.')
    parser_modifygenomelist.add_argument('--name', dest = 'name',
                                        help='Modify the name of the list to this.')
    parser_modifygenomelist.add_argument('--public', dest = 'public', type=bool,
                                        help='Change whether the list is private or public.')
    parser_modifygenomelist.set_defaults(func=ModifyGenomeList)

# --------- Clone A Genome List

    parser_clonegenomelist = subparsers.add_parser('CloneGenomeList',
                                        help='Create a genome list from a list of accessions')
    parser_clonegenomelist.add_argument('--list_id', dest = 'list_id', type=int,
                                       required=True, help='File containing list of accessions')
    parser_clonegenomelist.add_argument('--name', dest = 'name',
                                       required=True, help='Name of the genome list')
    parser_clonegenomelist.add_argument('--description', dest = 'description',
                                       required=True, help='Brief description of the genome list')
    parser_clonegenomelist.add_argument('--public', dest = 'public', default=False,
                                       action='store_true', help='Make the list visible to all users.')
    parser_clonegenomelist.set_defaults(func=CloneGenomeList)


# --------- Delete A Genome List

    parser_deletegenomelist = subparsers.add_parser('DeleteGenomeList',
                                        help='Create a genome list from a list of accessions')
    parser_deletegenomelist.add_argument('--list_id', dest = 'list_id', type=int,
                                       required=True, help='ID of the genome list to delete')
    parser_deletegenomelist.add_argument('--force', dest = 'force', action='store_true',
                                        help='Do not prompt for confirmation of deletion')
    parser_deletegenomelist.set_defaults(func=DeleteGenomeList)

# -------- Show All Genome Lists

    parser_showallgenomelists = subparsers.add_parser('ShowAllGenomeLists',
                                        help='Create a genome list from a list of accessions')
    parser_showallgenomelists.add_argument('--owned', dest = 'self_owned',  default=False,
                                        action='store_true', help='Only show genome lists owned by you.')
    parser_showallgenomelists.set_defaults(func=ShowAllGenomeLists)

# -------- Generate Tree Data

    parser_createtreedata = subparsers.add_parser('CreateTreeData',
                                        help='Generate data to create genome tree')
    parser_createtreedata.add_argument('--core_lists', dest = 'core_lists', choices=('private', 'public', 'both'),
                                        help='Include the genomes from one or all of the ACE core genome lists in the output files.')
    parser_createtreedata.add_argument('--list_ids', dest = 'list_ids',
                                        help='Create genome tree data from these lists (comma separated).')
    parser_createtreedata.add_argument('--set_id', dest = 'marker_set_id',
                                        required=True, help='Use this marker set for the genome tree.')
    parser_createtreedata.add_argument('--tree_ids', dest = 'tree_ids',
                                        help='Add these tree_ids to the output, useful for including outgroups (comma separated).')
    parser_createtreedata.add_argument('--output', dest = 'out_dir',
                                        required=True, help='Directory to output the files')
    parser_createtreedata.add_argument('--profile', dest = 'profile',
                                        help='Marker profile to use (default: %s)' % (profiles.ReturnDefaultProfileName(),))
    parser_createtreedata.add_argument('--profile_args', dest = 'profile_args',
                                        help='Arguments to provide to the profile')
    parser_createtreedata.set_defaults(func=CreateTreeData)

# -------- Marker management subparsers

    parser_recalculatemarkers = subparsers.add_parser('RecalculateMarkers',
                                help='Recalculate markers')
    parser_recalculatemarkers.add_argument('--tree_ids', dest = 'tree_ids',
                                         help='List of Tree IDs (comma separated)')
    parser_recalculatemarkers.add_argument('--filename', dest = 'listfile',
                                         help='File containing list of Tree IDs (newline separated)')
    parser_recalculatemarkers.set_defaults(func=RecalculateMarkers)

    parser_recalculateallmarkers = subparsers.add_parser('RecalculateAllMarkers',
                                help='Recalculate all the markers')

    parser_recalculateallmarkers.set_defaults(func=RecalculateAllMarkers)

#--------- Metadata managements

    parser_addcustommetadata = subparsers.add_parser('AddCustomMetadata',
                                  help='Add custom metadata to the database')
    parser_addcustommetadata.add_argument('--xml_path', dest = 'xml_path',
                                        required=True, help='XML path of metadata to be added (e.g "data/custom/metadatafield)')
    parser_addcustommetadata.add_argument('--metadata_file', dest = 'metadata_file',
                                        required=True, help='File (tab separated) containing tree ids and data to be added at specified XML path')
    parser_addcustommetadata.set_defaults(func=AddCustomMetadata)


    parser_updatetaxonomies = subparsers.add_parser('UpdateTaxonomies',
                                        help='Update the internal taxonomies')
    parser_updatetaxonomies.add_argument('--taxonomy_file', dest = 'taxonomy_file',
                                        required=True, help='File containing tree ids and taxonomies (tab separated)')
    parser_updatetaxonomies.set_defaults(func=UpdateTaxonomies)

#--------- Metadata managements - ModifyCoreLists

    parser_modifycorelist = subparsers.add_parser('ModifyCoreLists',
                                help='Add/remove genomes the private/public core list.')
    parser_modifycorelist.add_argument('--tree_ids', dest = 'tree_ids',
                                         required=True,  help='List of Tree IDs (comma separated)')
    parser_modifycorelist.add_argument('--operation', dest = 'operation', choices=('private','public','remove'),
                                         required=True,  help='Operation to perform')
    parser_modifycorelist.set_defaults(func=ModifyCoreLists)

#--------- Marker Management

    parser_addmarkers = subparsers.add_parser('AddMarkers',
                                 help='Add in one or many marker HMMs into the database')
    parser_addmarkers.add_argument('--database_name', dest='dbname', required=True,
                                help='Name of the database that the markers belong to')
    parser_addmarkers.add_argument('--file', dest='file', required=True,
                                help='File containing the HMM model(s) for the marker(s)')
    parser_addmarkers.add_argument('--use_existing', dest='use_existing', action='store_true', default=False,
                                help='If copies of this marker already exist in the database, use the copy in the database. These markers will still be added to the specified marker set (if specified)')
    mutex_group = parser_addmarkers.add_mutually_exclusive_group(required=False)
    mutex_group.add_argument('--modify_set', dest = 'marker_set_id',
                                    help='Modify a marker set with the \
                                    specified id and add all markers to it.')
    mutex_group.add_argument('--create_set', dest = 'marker_set_name',
                                    help='Create a marker set with the specified name and add these markers to it.')
    parser_addmarkers.set_defaults(func=AddMarkers)

    parser_showallmarkerdatabases = subparsers.add_parser('ShowAllMarkerDatabases',
                                                          help='Shows all the possible databases markers can com from')
    parser_showallmarkerdatabases.set_defaults(func=ModifyCoreLists)

#--------- Marker Management

    parser_deletemarkers = subparsers.add_parser('DeleteMarkers',
                                    help='Delete markers from the database')
    parser_deletemarkers.add_argument('--marker_ids', dest = 'marker_ids',
                                    required=True, help='List of Marker IDs (comma separated)')
    parser_deletemarkers.set_defaults(func=DeleteMarkers)

# -------- Show All Genome Lists

    parser_showallmarkersets = subparsers.add_parser('ShowAllMarkerSets',
                                        help='Shows the details of a Marker Set')
    parser_showallmarkersets.set_defaults(func=ShowAllMarkerSets)
    """

    args = parser.parse_args()

    # Special parser checks

    if (args.category_parser_name == 'trees' and args.tree_subparser_name == 'create'):
        if (args.genome_batchfile is None and args.genome_ids is None and args.genome_list_ids is None):
            parser_tree_create.error('Need to specify at least one of --genome_batchfile, --genome_ids or --genome_list_ids')
        if (args.marker_batchfile is None and args.marker_ids is None and args.marker_set_ids is None ):
            parser_tree_create.error('Need to specify at least one of --marker_batchfile, --marker_ids or --marker_set_ids')

    if (args.category_parser_name == 'genomes' and args.genome_subparser_name == 'view'):
        if (args.batchfile is not None or args.id_list is not None):
            if args.view_all:
                parser_genome_view.error('argument --all must be used by itself')
        elif not args.view_all:
            parser_genome_view.error('need to specify at least one of --all, --batchfile or --genome_ids')
            
    if (args.category_parser_name == 'genomes' and args.genome_subparser_name == 'delete'):
        if (args.batchfile is None and args.id_list is None):
                parser_genome_delete.error('need to specify at least one of --batchfile or --genome_ids')
                
    if (args.category_parser_name == 'markers' and args.marker_subparser_name == 'view'):
        if (args.batchfile is not None or args.id_list is not None):
            if args.view_all:
                parser_marker_view.error('argument --all must be used by itself')
        elif not args.view_all:
            parser_marker_view.error('need to specify at least one of --all, --batchfile or --marker_ids')

    # Initialise the backend
    db = GenomeDatabase.GenomeDatabase()
    db.conn.MakePostgresConnection()

    if args.debug:
        db.SetDebugMode(True)

    # Login
    if args.login_as_root:
        user = db.RootLogin(GetLinuxUsername())
    else:
        user = db.UserLogin(GetLinuxUsername())

    if not user:
        ErrorReport("Database login failed. The following error(s) were reported:\n")
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
        ErrorReport("Database action failed. The following error(s) were reported:\n")
        DumpDBErrors(db)
        sys.exit(-1)

