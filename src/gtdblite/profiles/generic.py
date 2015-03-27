import os
import sys
import psycopg2 as pg
import common

valid_configs = [('individual', type(None), "Create individual FASTA files for each marker instead of a concatenated alignment."),
                 ('taxonomy', type(''), "Filter for organisms with this taxonomy (internal genome tree taxonomy)"),
                 ('reject_missing_taxonomy', type(None), "Reject any genomes that have not got an assigned taxonomy.")]

def GetValidConfigOptions():
    return valid_configs
    
def MakeTreeData(GenomeDatabase, marker_ids, genome_ids, directory, prefix=None, config_dict=None):
    if not os.path.isdir(directory):
        GenomeDatabase.ReportError("Directory doesn't exist: " + directory)
        return None
    
    if not common.CheckPassedConfigsAgainstKnownConfigs(config_dict, GetValidConfigOptions()):
        return None
    
    cur = GenomeDatabase.conn.cursor()
      
    # Get the total number of markers in this marker set
    total_marker_count = len(marker_ids)
    
    # For each genome, get the number of markers in this marker set that genome contains    
    genome_gene_counts = GenomeDatabase.GetAlignedMarkersCountForGenomes(genome_ids, marker_ids)
    
    # For all of the markers, get the expected marker size.
    cur.execute("SELECT id, id_in_database, size " +
                "FROM markers " +
                "WHERE id in %s " 
                "ORDER by markers.id", (tuple(marker_ids),))
    
    chosen_markers = dict()
    for marker_id, id_in_database, size in cur:
        chosen_markers[marker_id] = {'id_in_database': id_in_database, 'size': size}
    
    individual_marker_fasta = dict()
    genome_info = dict()
    
    if prefix is None:
        prefix = "genome_tree_data"

    gg_fh = open(os.path.join(directory, prefix + "_concatenated.greengenes"), 'wb')
    fasta_concat_fh = open(os.path.join(directory, prefix + "_concatenated.faa"), 'wb')
    
    # Run through each of the genomes and make the magic happen.
    for genome_id in genome_ids:
        
        cur.execute("SELECT external_id_prefix || '_' || id_at_source as external_id, genomes.name, owned_by_root, users.username "+
                    "FROM genomes "+
                    "LEFT OUTER JOIN users ON genomes.owner_id = users.id " +
                    "LEFT OUTER JOIN genome_sources ON genomes.genome_source_id = genome_sources.id " +
                    "WHERE genomes.id = %s", (genome_id,))
        
        result = cur.fetchone()
        if not result:
            continue
        (external_id, name, owned_by_root, owner) = result
        
        
        # Populate genome info
        genome_info['markers']  = dict()
        genome_info['name']  = name
        genome_info['external_id']  = external_id
        genome_info['owner']  = ('root' if owned_by_root else owner)
            
        cur.execute("SELECT aligned_markers.marker_id, sequence " +
                    "FROM aligned_markers "+
                    "WHERE genome_id = %s " +
                    "AND sequence is NOT NULL "+
                    "AND marker_id in %s ", (genome_id, tuple(marker_ids,)))
        if (cur.rowcount == 0):
            sys.stderr.write("WARNING: Genome %s has no markers for this marker set in the database and will be missing from the output files.\n" % external_id)
            sys.stderr.flush()
            continue
        for marker_id, sequence in cur:
            genome_info['markers'][marker_id] = sequence
            
        aligned_seq = '';
        for marker_id in chosen_markers.keys():
            if marker_id in genome_info['markers']:
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
    
        gg_list = ["BEGIN",
                    "db_name=%s" % genome_info['external_id'],
                    "organism=%s" % genome_info['name'],
                    "prokMSA_id=%s" % genome_info['external_id'],
                    "owner=%s" % genome_info['owner'],
                    "remark=%iof%i" % (genome_gene_counts[genome_id], total_marker_count),
                    "warning=",
                    "aligned_seq=%s" % (aligned_seq),
                    "END"]
        
        gg_outstr = "\n".join(gg_list) + "\n\n";
        
        fasta_concat_fh.write(fasta_outstr)
        gg_fh.write(gg_outstr)
    
    gg_fh.close()
    fasta_concat_fh.close()


    if "individual" in config_dict:
        for marker_id in chosen_markers.keys():
            fasta_individual_fh = open(os.path.join(directory, prefix + "_" + str(marker_id) + "_" + 
                                                               chosen_markers[marker_id]['id_in_database'] + ".faa"),
                                           'wb')
            fasta_individual_fh.write(''.join(individual_marker_fasta[marker_id]))
            fasta_individual_fh.close()
    
    
    return True