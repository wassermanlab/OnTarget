from os.path import dirname, abspath, join
from os import listdir
from json import loads

def get_sample_names_and_uids_from_organs(args):
    try:
        organs = args['organs']
        optional_args = _get_optional_args(args)
        organ_json_file_paths = _get_organ_json_file_paths(organs)
        matching_sample_names_and_uids = _get_list_of_matching_sample_names_and_uids(organ_json_file_paths, optional_args)
    except:
        return None
    
    return {'result': matching_sample_names_and_uids}


def _get_list_of_matching_sample_names_and_uids(organ_json_file_paths, optional_args):
    matching_sample_names_and_uids = []

    for organ_json_file_path in organ_json_file_paths:
        organ_sample_dict = _json_file_path_to_dict(organ_json_file_path)

        for samples in organ_sample_dict.values():
            for sample in samples:
                if _sample_matches_optional_args(sample, optional_args):
                    matching_sample_names_and_uids.append({'name': sample['name'], 'uid': sample['uid']})
    
    return matching_sample_names_and_uids


def _json_file_path_to_dict(json_file_path):
    organ_json_string = open(json_file_path).read()
    organ_sample_dict = loads(organ_json_string)

    return organ_sample_dict


def _sample_matches_optional_args(sample, optional_args):
    for key, value in optional_args.items():
        if value is not None:
            if sample[key] != value:
                return False

    return True    


def _get_organ_json_file_paths(organs):
    unique_organs = set(organs)
    organ_json_file_paths = []

    root_path = abspath(dirname(dirname(abspath(__file__))))
    samples_path = join(root_path, 'samples')
    sample_file_names = listdir(samples_path)

    for file_name in sample_file_names:
        if file_name.endswith('.json'):
            if file_name[:len(file_name) - 5] in unique_organs:
                organ_json_file_path = join(samples_path, file_name)
                organ_json_file_paths.append(organ_json_file_path)
    
    if len(organ_json_file_paths) != len(unique_organs):
        raise Exception("Unable to locate sample json files for all organs")
    
    return organ_json_file_paths


def _get_optional_args(args):
    optional_args = {
        'X': None,
        'Y': None,
        'cancer': None,
        'cell_line': None,
        'treatment': None
    }

    for key in optional_args.keys():
        if hasattr(args, key):
            optional_args[key] = args[key]

    return optional_args