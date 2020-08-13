import traceback

import requests
from fuzzywuzzy import process, fuzz


def name2sample(remote_host, database, name=None):
    URL = remote_host + '/api/v1/' + database + '/samples'
    # name = request.args.get('name')
    results = {}
    try:
        data = requests.get(url=URL).json()
        print(data)
        list = []
        while "next" in data:
            list.extend(data['results'])
            data = requests.get(url=data['next']).json()
        list.extend(data['results'])
        samples = {}
        for item in list:
            samples.setdefault(item['name'], [])
            samples[item["name"]].append(item)
        if name is not None:
            fuzzy = process.extract(name, samples.keys(), scorer=fuzz.token_set_ratio,
                                    limit=None)
            formatted_results = []
            for name, score in fuzzy:
                for sample in samples[name]:
                    sample.setdefault("relevance", score)
                    formatted_results.append(sample)

            results['results'] = formatted_results
        else:
            results['results']=list
    except:
        traceback.print_exc()
        results["Error"] = "Could not get samples"
    return results
