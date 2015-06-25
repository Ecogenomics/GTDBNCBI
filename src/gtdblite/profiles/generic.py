import os
import sys
import psycopg2 as pg
import common

valid_configs = [('individual', type(None), "Create individual FASTA files for each marker instead of a concatenated alignment."),
                 ('checkm_contamination_threshold', float, "Only include genomes with CheckM contamination below this"),
                 ('checkm_completeness_threshold', float, "Only include genomes with CheckM completeness above this"),
                 ('include_multihits', type(None), "By default, aligned genes that have multiple vetted alignment are excluded from the alignment. \
                  This flag tells the profile to use the best hit of the multiple hits (may create chimeric sequences)"),
                 ('guaranteed_genome_ids', str, "Comma separated list of genome IDs that will not be filtered and a guaranteed to be placed in the tree."),
                 ('guaranteed_genome_list_ids', str, "Comma separated list of genome list IDs whose genomes will not be filtered and a guaranteed to be placed in the tree."),] 

default_checkm_contamination = 10
default_checkm_completeness = 50

def GetValidConfigOptions():
    return valid_configs
    
def MakeTreeData(GenomeDatabase, marker_ids, genome_ids, directory, prefix=None, config_dict=None):
    if not os.path.isdir(directory):
        GenomeDatabase.ReportError("Directory doesn't exist: " + directory)
        return None
    
    if not common.CheckPassedConfigsAgainstKnownConfigs(GenomeDatabase, config_dict, GetValidConfigOptions()):
        return None
    
    if 'checkm_contamination_threshold' not in config_dict:
        GenomeDatabase.ReportWarning("'checkm_contamination_threshold' flag not supplied to tree building profile 'generic'. Using default cutoff of %f." % default_checkm_contamination)
        config_dict['checkm_contamination_threshold'] = default_checkm_contamination
    
    if 'checkm_completeness_threshold' not in config_dict:
        GenomeDatabase.ReportWarning("'checkm_completeness_threshold' flag not supplied to tree building profile 'generic'. Using default cutoff of %f." % default_checkm_completeness)
        config_dict['checkm_completeness_threshold'] = default_checkm_completeness
    
    cur = GenomeDatabase.conn.cursor()

    # For all of the markers, get the expected marker size.
    if GenomeDatabase.debugMode:
        sys.stderr.write("Getting expected marker sizes...\n")
        sys.stderr.flush()    
    
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
    
    individual_marker_fasta = dict()
    genome_info = dict()
    
    if prefix is None:
        prefix = "genome_tree_data"

    gg_fh = open(os.path.join(directory, prefix + "_concatenated.arbtxt"), 'wb')
    fasta_concat_fh = open(os.path.join(directory, prefix + "_concatenated.faa"), 'wb')
    marker_info_fh = open(os.path.join(directory, prefix + "_concatenated_markers.info"), 'wb')
    multi_hits_fh = open(os.path.join(directory, prefix + "_multi_hits.info"), 'wb')
    
    
    # Output the marker info and multi_hit info
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
    
    
    # Find genomes that are in the guaranteed list
    guaranteed_genomes = set()
    if 'guaranteed_genome_ids' in config_dict:
        guaranteed_genomes.update(GenomeDatabase.ExternalGenomeIdsToGenomeIds(config_dict['guaranteed_genome_ids'].split(",")))
        
    if 'guaranteed_genome_list_ids' in config_dict:
        guaranteed_genomes.update(GenomeDatabase.GetGenomeIdListFromGenomeListIds(config_dict['guaranteed_genome_list_ids'].split(",")))
    
    if GenomeDatabase.debugMode:
        sys.stderr.write("Running checkM filters...\n")
        sys.stderr.flush()      
    
    # Find genomes that fail checkm cutoffs
    cur.execute("SELECT external_id_prefix || '_' || id_at_source as external_id, genomes.id, checkm_completeness, checkm_contamination "+
                 "FROM genomes "+
                 "LEFT OUTER JOIN genome_sources ON genomes.genome_source_id = genome_sources.id " +
                 "WHERE genomes.id in %s "+
                 "AND (checkm_completeness < %s "+
                      "OR checkm_contamination > %s)",
                (tuple(genome_ids), config_dict['checkm_completeness_threshold'], config_dict['checkm_contamination_threshold']))
    
    filtered_genomes = set()
    for (external_id, genome_id, checkm_completeness, checkm_contamination) in cur:
        filtered_genomes.add(genome_id)
    
    
    if GenomeDatabase.debugMode:
        sys.stderr.write("Outputting aligned markers...\n")
        sys.stderr.flush()  
    
    # Select all the genomes - could be lots of genome id, create a temp table
    temp_table_name = GenomeDatabase.GenerateTempTableName()

    cur.execute("CREATE TEMP TABLE %s (id integer)" % (temp_table_name,) )
    query = "INSERT INTO {0} (id) VALUES (%s)".format(temp_table_name)
    cur.executemany(query, [(x,) for x in genome_ids])
    
    query = ("SELECT external_id_prefix || '_' || id_at_source as external_id, genomes.id, genomes.name, checkm_completeness, checkm_contamination, owned_by_root, users.username " +
             "FROM genomes "+
             "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
             "LEFT OUTER JOIN genome_sources ON genomes.genome_source_id = genome_sources.id " +
             "WHERE genomes.id in (SELECT id from {0})".format(temp_table_name))
    
    cur.execute(query) 

    genome_details = cur.fetchall()
    
    # Run through each of the genomes and make the magic happen.
    for (external_id, genome_id, name, checkm_completeness, checkm_contamination, owned_by_root, owner) in genome_details:
        
        if genome_id in filtered_genomes:
            
            if genome_id not in guaranteed_genomes:
                sys.stderr.write("WARNING: Genome id %s excluded. Failed checkM completeness/contamination cutoffs. This genome's values: (%f, %f)\n" %
                                 (external_id, checkm_completeness, checkm_contamination))
                sys.stderr.flush()
                continue

        # Populate genome info
        genome_info['markers']  = dict()
        genome_info['multiple_hits']  = dict()
        genome_info['name']  = name
        genome_info['external_id']  = external_id
        genome_info['owner']  = ('root' if owned_by_root else owner)
            
        cur.execute("SELECT aligned_markers.marker_id, sequence, multiple_hits " +
                    "FROM aligned_markers "+
                    "WHERE genome_id = %s " +
                    "AND sequence is NOT NULL "+
                    "AND marker_id in %s ", (genome_id, tuple(marker_ids,)))
        
        if (cur.rowcount == 0):
            sys.stderr.write("WARNING: Genome %s has no markers for this marker set in the database and will be missing from the output files.\n" % external_id)
            sys.stderr.flush()
            continue
        
        for marker_id, sequence, multiple_hits in cur:
            genome_info['markers'][marker_id] = sequence
            genome_info['multiple_hits'][marker_id] = multiple_hits
        
        aligned_seq = '';
        multi_hits_details = [genome_info['external_id']]
        for marker_id in chosen_markers_order:
            multiple_hits = genome_info['multiple_hits'][marker_id]
            multi_hits_details.append(
                'Multiple' if multiple_hits else 'Single'
            )
            if (marker_id in genome_info['markers']) and (not multiple_hits or 'include_multihits' in config_dict):
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
                    "multiple_homologs=%i/%i" % (sum([1 if x == "Multiple" else 0 for x in multi_hits_details]), len(chosen_markers_order)),# Lazy, lazy. Fix this.
                    "warning=",
                    "aligned_seq=%s" % (aligned_seq),
                    "END"]
        
        gg_outstr = "\n".join(gg_list) + "\n\n";
        
        fasta_concat_fh.write(fasta_outstr)
        gg_fh.write(gg_outstr)
        multi_hits_fh.write(multi_hits_outstr)
    
    gg_fh.close()
    fasta_concat_fh.close()
    multi_hits_fh.close()
    
    if "individual" in config_dict:
        for marker_id in chosen_markers.keys():
            fasta_individual_fh = open(os.path.join(directory, prefix + "_" + str(marker_id) + "_" + 
                                                               chosen_markers[marker_id]['id_in_database'] + ".faa"),
                                           'wb')
            fasta_individual_fh.write(''.join(individual_marker_fasta[marker_id]))
            fasta_individual_fh.close()
    
    
    return True
