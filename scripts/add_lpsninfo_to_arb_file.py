#!/usr/bin/env python
import argparse
import os
import sys

def warn(msg):
    sys.stderr.write(str(msg))
    sys.stderr.flush()

def warnln(msg):
    sys.stderr.write(str(msg) + "\n")
    sys.stderr.flush()

def add_tax_to_arbfile(lpsn_info_file, arb_file):
    
    map_genus_to_info = {}
    map_species_to_info = {}
    
    info_fh = open(lpsn_info_file, "rb")
    for line in info_fh:
        splitline = line.rstrip().split("\t")
        full_name = splitline[2] 
        if splitline[0] == 'genus':
            map_genus_to_info[full_name] = {'info': splitline[3]};
        elif splitline[0] == 'species':
            
            splitname = full_name.split()
            genus = splitname[0]
            type_species = (True if splitline[1] == "1" else False)
            strains = splitline[4:]
            map_species_to_info[full_name] = {'info': splitline[3], 'strains': strains, 'type_species': type_species, 'genus': genus}
    info_fh.close()
    
    map_previously_added_species = {}
    
    current_id = None
    species_name = None
    genus_name = None
    
    arb_fh = open(arb_file, "rb")
    for line in arb_fh:
        line = line.rstrip()
        splitline = line.split("=")
        if line[0:3] == "END":
            current_id = None
        elif splitline[0] == "db_name":
            current_id = splitline[1]
        elif splitline[0] == "organism":
            organism_name = splitline[1]
            species_name = None
            
            organism_split_name = organism_name.split()
            genus_name = organism_split_name[0]
            
            # Do some checks
            
            matches = []
            for this_known_species in map_species_to_info.keys():
                if organism_name.find(this_known_species) != -1: 
                    matches.append(this_known_species)

            if not matches:
                organism_name = None
                
            if len(matches) > 1:
                
                longest_match = ''
                for match in matches:
                    if len(match) > len(longest_match):
                        longest_match = match
                warnln("Multiple species hits for %s. Matches: %s. Will use %s." % (organism_name, str(matches), longest_match))
                matches = [longest_match]
                
            
            if organism_name:
            
                species_name = matches[0]
                
                species_name_start = organism_name.find(species_name)
                
                matches = []
                for strain_name in map_species_to_info[species_name]['strains']:
                    if len(strain_name) < 2:
                        warnln("Strain name %s from species %s is too short. May result in ambiguous hits. Ignoring." % (strain_name, species_name))
                        continue
                    if organism_name.find(strain_name, species_name_start + len(species_name)) != -1: 
                        matches.append(strain_name)
                        
                if not matches:
                    species_name = None

                if len(matches) > 1:
                    
                    longest_match = ''
                    for match in matches:
                        if len(match) > len(longest_match):
                            longest_match = match
                    warnln("Multiple strain hits for %s. Matches: %s. Will use %s." % (organism_name, str(matches), longest_match))
                    matches = [longest_match]

                if species_name:
                    strain_hit_string = matches[0]
                    
                    this_strain_info = {'organism_name': organism_name, 'strain': strain_hit_string, 'db_id': current_id}
    
                    if species_name in map_previously_added_species:
                        warnln("Multiple type strains found for %s.  This strain info: %s. Previous strain info: %s." % (species_name, this_strain_info, str(map_previously_added_species[species_name])))
                    
                    map_previously_added_species[species_name] = this_strain_info

                    
        elif (splitline[0] == "warning"):
            if genus_name is not None:
                try:
                    print 'genus_info=(%s)%s' % (genus_name, map_genus_to_info[genus_name]['info'])
                except KeyError:
                    pass
                
            if (current_id is not None) and (species_name is not None):
                information = []
                if map_species_to_info[species_name]['type_species']:
                    genus = map_species_to_info[species_name]['genus']
                    information.append('type_species (Info: %s)' % (map_genus_to_info[genus]['info']))
                
                information.append('type_strain (Match String: "%s"; Info: %s)' % (map_previously_added_species[species_name]['strain'], map_species_to_info[species_name]['info']))
                
                print 'type_material=%s' % "; ".join(information)
            
            
        print line
    arb_fh.close()
 
if __name__ == '__main__':

    # create the top-level parser
    parser = argparse.ArgumentParser(prog='add_lpsninfo_to_arb_file.py')
    parser.add_argument('-l', dest='lpsn_info_file', required=True,
                        help='File containing lpsn info')
    parser.add_argument('-a', dest='arbfile', required=True, 
                        help='Arb txt file')
    args = parser.parse_args()
    
    add_tax_to_arbfile(args.lpsn_info_file, args.arbfile) 