#!/usr/bin/env python
import argparse
import sys
from gtdblite import GenomeDatabase 
import os
import pwd

def GetLinuxUsername():
    return pwd.getpwuid(os.getuid())[0]

def ErrorReport(msg):
    sys.stderr.write(msg)
    sys.stderr.flush()

def AddManyFastaGenomes(db, args):
    return db.AddManyFastaGenomes(args.batchfile, args.checkm_file, args.force)

def AddMarkers(db, args):
    return db.AddMarkers(args.batchfile, args.force)

def CreateTreeData(db, args):    
        
        list_ids = []
        if args.list_ids:
            list_ids = args.list_ids.split(",")
            
        genome_id_set = set()
        if args.tree_ids:
            extra_ids = [db.GetGenomeId(x) for x in args.tree_ids.split(",")]
            genome_id_set = genome_id_set.union(set(extra_ids))
        for list_id in list_ids:
            temp_genome_list = db.GetGenomeIdListFromGenomeListId(list_id)
            if temp_genome_list:
                genome_id_set = genome_id_set.union(set(temp_genome_list))

        profile_config_dict = dict()
        if args.profile_args:
            profile_args = args.profile_args.split(',')
            for profile_arg in profile_args:
                key_value_pair = profile_arg.split('=')
                try:
                    profile_config_dict[key_value_pair[0]] = key_value_pair[1]
                except IndexError:
                    profile_config_dict[key_value_pair[0]] = None
        
        if (len(genome_id_set) > 0) or (len(core_lists) != 0):
            db.MakeTreeData(args.marker_set_id,  list(genome_id_set), args.profile, args.out_dir, config_dict=profile_config_dict)

if __name__ == '__main__':
    
    # create the top-level parser
    parser = argparse.ArgumentParser(prog='gtdblite.py')
    parser.add_argument('-r', dest='login_as_root', action='store_true',
                        help='Login as the root user'),
    parser.add_argument('-f', dest='force', action='store_true', 
                        help='Force the action (required to override warnings for certain actions)'),
    parser.add_argument('--dev', dest='dev', action='store_true',
                        help='Connect to the developer database')
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Run in debug mode')
    
    category_parser = parser.add_subparsers(help='Category Command Help', dest='category_parser_name')

    user_category_parser = category_parser.add_parser('user', help='View the user sub-commands')
    user_category_subparser = user_category_parser.add_subparsers(help='User command help', dest='user_subparser_name')
    
    genome_category_parser = category_parser.add_parser('genome', help='View the genome sub-commands')
    genome_category_subparser = genome_category_parser.add_subparsers(help='Genome command help', dest='genome_subparser_name')
    
    marker_category_parser = category_parser.add_parser('marker', help='View the marker commands')
    marker_category_subparser = marker_category_parser.add_subparsers(help='Marker command help', dest='marker_subparser_name')
    
    tree_category_parser = category_parser.add_parser('tree', help='View the tree commands')
    tree_category_subparser = tree_category_parser.add_subparsers(help='Tree command help', dest='tree_subparser_name')
    
    profile_category_parser = category_parser.add_parser('profile', help='View the profile commands')
    profile_category_subparser = profile_category_parser.add_subparsers(help='Profile command help', dest='profile_subparser_name')
# -------- Genome Management subparsers

    parser_addmanyfastagenomes = genome_category_subparser.add_parser('add',
                                    help='Add one or many genomes to the tree.')
    parser_addmanyfastagenomes.add_argument('--batchfile', dest = 'batchfile',
                                    required=True, help='Batchfile describing the genomes - one genome per line, tab separated in 3-5 columns (filename, name, desc, [source], [id_at_source])')
    parser_addmanyfastagenomes.add_argument('--checkm_result', dest = 'checkm_file',
                                    required=True, help='Provide a checkM results file. MUST BE A TAB TABLE! e.g. "checkm taxonomy_wf -f CHECKM_FILE --tab_table domain Bacteria bins/ output"')
    mutex_group = parser_addmanyfastagenomes.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--modify_list', dest = 'genome_list_id',
                                    help='Modify a genome list with the \
                                    specified id and add all batchfile genomes into it.')
    mutex_group.add_argument('--create_list', dest = 'genome_list_name',
                                    help='Create a genome list with the specified name and add all batchfile genomes into it.')
    mutex_group.add_argument('--no_list', dest = 'no_genome_list', action="store_true",
                                    help="Don't add these genomes to a list.")
    
    parser_addmanyfastagenomes.set_defaults(func=AddManyFastaGenomes)

#--------- Marker Management Subparsers
    
    parser_addmarkers = marker_category_subparser.add_parser('add', 
                                 help='Add in one or many marker HMMs into the database')
    parser_addmarkers.add_argument('--batchfile', dest='batchfile', required=True,
                                help='Batchfile describing the markers - one HMM file per line (one model per file), tab separated in 3-5 columns (filename, name, desc, [database], [database_specific_id]')
    mutex_group = parser_addmarkers.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--modify_set', dest = 'marker_set_id',
                                    help='Modify a marker set with the \
                                    specified id and add all markers to it.')
    mutex_group.add_argument('--create_set', dest = 'marker_set_name',
                                    help='Create a marker set with the specified name and add these markers to it.')
    mutex_group.add_argument('--no_set', dest = 'no_marker_set', action="store_true",
                                    help="Don't add these markers to a marker set.")
    parser_addmarkers.set_defaults(func=AddMarkers)

# -------- Generate Tree Data

    ## FROM HERE!!!

    parser_createtreedata = tree_category_subparser.add_parser('create',
                                        help='Generate data to create genome tree')
    parser_createtreedata.add_argument('--list_id_file', dest = 'list_id_file',
                                        help='Create genome tree data from these lists (comma separated).')
    parser_createtreedata.add_argument('--set_id', dest = 'marker_set_id',
                                        required=True, help='Use this marker set for the genome tree.')
    parser_createtreedata.add_argument('--output', dest = 'out_dir',
                                        required=True, help='Directory to output the files')
    parser_createtreedata.add_argument('--profile', dest = 'profile',
                                        help='Tree creation profile to use (default: %s)' % (profiles.ReturnDefaultProfileName(),))
    parser_createtreedata.add_argument('--profile_args', dest = 'profile_args',
                                        help='Arguments to provide to the profile')
    parser_createtreedata.set_defaults(func=CreateTreeData)

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
    
    # Initialise the backend
    db = GenomeDatabase.GenomeDatabase()
    db.conn.MakePostgresConnection(args.dev)

    if args.debug:
        db.SetDebugMode(True)
    
    # Login
    if args.login_as_root:
        user = db.RootLogin(GetLinuxUsername())    
    else:
        user = db.UserLogin(GetLinuxUsername())
        
    if not user:
        ErrorReport("Database login failed. The following error was reported:\n" +
                    "\t" + db.lastErrorMessage)
        sys.exit(-1)

    result = args.func(db, args)
    
    if not result:
        ErrorReport("Database action failed. The following error was reported:\n" +
            "\t" + db.lastErrorMessage)
    
    
    
    
    
    