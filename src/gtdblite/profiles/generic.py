import os
import sys
import psycopg2 as pg
import common

valid_configs = [('individual', type(None), "Create individual FASTA files for each marker instead of a concatenated alignment."),
                 ('taxonomy', type(''), "Filter for organisms with this taxonomy (internal genome tree taxonomy)"),
                 ('reject_missing_taxonomy', type(None), "Reject any genomes that have not got an assigned taxonomy.")]

def GetValidConfigOptions():
    return valid_configs
    
def MakeTreeData(GenomeDatabase, marker_set_id, list_of_genome_ids, directory, prefix=None, config_dict=None):
    if not os.path.isdir(directory):
        GenomeDatabase.ReportError("Directory doesn't exist: " + directory)
        return None
    
    if not common.CheckPassedConfigsAgainstKnownConfigs(config_dict, GetValidConfigOptions()):
        return None
    
    cur = GenomeDatabase.conn.cursor()
      
    # Get the total number of markers in this marker set
    total_marker_count = len(GenomeDatabase.GetMarkerIdListFromMarkerSetId(marker_set_id))
    
    # For each genome, get the number of markers in this marker set that genome contains    
    genome_gene_counts = GenomeDatabase.GetAlignedMarkersCountForGenomeFromMarkerSetId(marker_set_id)
    
    # This filter needs to know about phylosift markers, calculate them if they haven't been calculated
    for genome_id in list_of_genome_ids:
        uncalculated_markers = GenomeDatabase.FindUncalculatedMarkersForGenomeId(genome_id, 
                                    GenomeDatabase.GetMarkerIdListFromMarkerSetId(1))
        if len(uncalculated_markers) > 0:
            if GenomeDatabase.debugMode:
                print "Calculating markers for ", genome_id
                print uncalculated_markers  
            GenomeDatabase.RecalculateMarkersForGenome(genome_id, uncalculated_markers)
        
    # For all of the markers, get the expected marker size.
    cur.execute("SELECT markers.id, markers.database_specific_id, size " +
                "FROM markers, marker_set_contents " +
                "WHERE set_id = %s " +
                "AND marker_id = markers.id "
                "ORDER by markers.id", (marker_set_id,))
    
    chosen_markers = dict()
    for marker_id, database_specific_id, size in cur:
        chosen_markers[marker_id] = {'database_specific_id': database_specific_id, 'size': size}
    
    individual_marker_fasta = dict()
    genome_info = dict()
    
    if prefix is None:
        prefix = "genome_tree_data"

    gg_fh = open(os.path.join(directory, prefix + "_concatenated.greengenes"), 'wb')
    fasta_concat_fh = open(os.path.join(directory, prefix + "_concatenated.faa"), 'wb')
    
    # Run through each of the genomes and make the magic happen.
    for genome_id in list_of_genome_ids:
        cur.execute("SELECT tree_id, name, XMLSERIALIZE(document metadata as text), username "+
                    "FROM genomes, users "+
                    "WHERE users.id = owner_id "+
                    "AND genomes.id = %s", (genome_id,))
        result = cur.fetchone()
        if not result:
            continue
        (tree_id, name, xmlstr, owner) = result
        
        # Populate genome info
        genome_info['markers']  = dict()
        genome_info['name']  = name
        genome_info['tree_id']  = tree_id
        genome_info['xmlstr']  = xmlstr
        genome_info['owner']  = owner
        
        # update XML metadata
        genome_info.update(common.GetInternalMetadataDictFromXMLString(xmlstr))
    
        # Test taxonomy
        taxonomy_success = True
        if not genome_info['internal_tax']:
            if 'reject_missing_taxonomy' in config_dict:
                taxonomy_success = False
        else:
            if ('taxonomy' in config_dict and config_dict['taxonomy']):
                if (genome_info['internal_tax'].find(config_dict['taxonomy']) < 0):
                   taxonomy_success = False
        if not taxonomy_success:
            if GenomeDatabase.debugMode:
                sys.stderr.write("%s (%s) not included. Filtered at taxonomy stage (taxonomy: %s)\n" % (genome_info['tree_id'],
                                                                                                      genome_info['name'],
                                                                                                      genome_info['internal_tax']))
                sys.stderr.flush()
            continue
            
        cur.execute("SELECT aligned_markers.marker_id, sequence " +
                    "FROM aligned_markers, marker_set_contents "+
                    "WHERE marker_set_contents.marker_id = aligned_markers.marker_id " +
                    "AND genome_id = %s " +
                    "AND sequence is NOT NULL "+
                    "AND set_id = %s ", (genome_id, marker_set_id))
        if (cur.rowcount == 0):
            sys.stderr.write("WARNING: Genome %s has no markers for this marker set in the database and will be missing from the output files.\n" % tree_id)
            sys.stderr.flush()
            continue
        for marker_id, sequence in cur:
            genome_info['markers'][marker_id] = sequence
            
        aligned_seq = '';
        for marker_id in chosen_markers.keys():
            if marker_id in genome_info['markers']:
                sequence = genome_info['markers'][marker_id]
                fasta_outstr = ">%s\n%s\n" % (genome_info['tree_id'],
                                              sequence)
                try:
                    individual_marker_fasta[marker_id].append(fasta_outstr)
                except KeyError:
                    individual_marker_fasta[marker_id] = [fasta_outstr]
            else:
                sequence = chosen_markers[marker_id]['size'] * '-'
            aligned_seq += sequence
        
        fasta_outstr = ">%s\n%s\n" % (genome_info['tree_id'],
                                      aligned_seq)
    
        gg_list = ["BEGIN",
                    "db_name=%s" % genome_info['tree_id'],
                    "organism=%s" % genome_info['name'],
                    "prokMSA_id=%s" % genome_info['tree_id'],
                    "owner=%s" % genome_info['owner'],
                    "genome_tree_tax_string=%s" % genome_info['internal_tax'],
                    "greengenes_tax_string=%s" % genome_info['gg_tax'],
                    "core_list_status=%s" % genome_info['core_list_status'],
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
                                                               chosen_markers[marker_id]['database_specific_id'] + ".faa"),
                                           'wb')
            fasta_individual_fh.write(''.join(individual_marker_fasta[marker_id]))
            fasta_individual_fh.close()
    
    
    return True