import xml.etree.ElementTree as ET
import sys

def GetInternalMetadataDictFromXMLString(xmlstr):
    root = ET.fromstring(xmlstr)
    extant = root.findall('internal/greengenes/dereplicated/best_blast/greengenes_tax')
    metadata = dict()
    metadata['gg_tax'] = ''
    if len(extant) != 0:
        metadata['gg_tax'] = extant[0].text
    
    extant = root.findall('internal/taxonomy')
    metadata['internal_tax'] = ''
    if len(extant) != 0:
        metadata['internal_tax'] = extant[0].text
    
    extant = root.findall('internal/core_list')
    metadata['core_list_status'] = ''
    if len(extant) != 0:
        metadata['core_list_status'] = extant[0].text
        
    return metadata

def ReportIncorrectParameter(db, param, expected_type, got_type, got_value):
    db.ReportError("Error: Incorrect parameter for config option: %s. Expected: %s. Got: %s (%s)" % (param, expected_type, got_type, got_value))
    sys.stderr.flush()
    
def CheckPassedConfigsAgainstKnownConfigs(db, passed_config_dict, valid_configs):
    valid_configs_type_dict = dict([(name, config_type) for (name, config_type, description) in valid_configs])
    for (config_key, value) in passed_config_dict.items():
        if config_key in valid_configs_type_dict:
            if valid_configs_type_dict[config_key] != type(value):
                if valid_configs_type_dict[config_key] == type(None):
                    db.ReportError("Unexpected value for valueless config option: %s. Got: %s" % (config_key, value))
                    return False
                elif valid_configs_type_dict[config_key] == type(1):
                    try:
                        int(value)
                    except ValueError:
                        ReportIncorrectParameter(db, config_key, valid_configs_type_dict[config_key], type(value), value)
                        return False
                elif valid_configs_type_dict[config_key] == type(0.1):
                    try:
                        float(value)
                    except ValueError:
                        ReportIncorrectParameter(db, config_key, valid_configs_type_dict[config_key], type(value), value)
                        return False
                else:
                    ReportIncorrectParameter(db, config_key, valid_configs_type_dict[config_key], type(value), value)
                    return False
        else:
            db.ReportWarning("Ignoring unknown config option (%s) for current profile." % (config_key))
            del passed_config_dict[config_key]
    return True