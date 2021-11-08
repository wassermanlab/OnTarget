from fuzzywuzzy import fuzz, process
import requests

from constants import GUD_URL

def get_sample_for_name(genome, name):
    try:
        list_of_samples = _get_list_of_samples_from_gud(genome)
        sample_name_to_samples_dict = _map_samples_to_their_names(list_of_samples)
        formatted_results = _get_formatted_results(name, sample_name_to_samples_dict)
    except:
        return None

    return {'result' : formatted_results}


def _get_list_of_samples_from_gud(genome):
    samples = []
    url = GUD_URL + genome + '/samples'

    while True:
        response = requests.get(url).json()
        samples.extend(response['results'])

        if 'next' in response:
            url = response['next']
        else:
            break
    
    return samples


def _map_samples_to_their_names(list_of_samples):
    name_to_samples_dict = {}

    for sample in list_of_samples:
        sample_name = sample['name']

        if sample_name in name_to_samples_dict:
            name_to_samples_dict[sample_name].append(sample)
        else:
            name_to_samples_dict[sample_name] = [sample]

    return name_to_samples_dict


def _get_formatted_results(name, sample_name_to_samples_dict):
    formatted_results = []

    results = process.extract(name, sample_name_to_samples_dict.keys(), 
        scorer=fuzz.token_set_ratio, limit=None)
    
    for name, score in results:
        for sample in sample_name_to_samples_dict[name]:
            sample.setdefault('relevance', score)
            formatted_results.append(sample)
    
    return formatted_results