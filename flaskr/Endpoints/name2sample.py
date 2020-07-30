import traceback

import requests
from fuzzywuzzy import process, fuzz


def name2sample(remote_host, database, name):
    URL = remote_host + '/api/v1/' + database + '/samples'
    # name = request.args.get('name')
    results = {}
    if name is not None:
        try:
            data = requests.get(url=URL).json()
            list = []
            while "next" in data:
                list.extend(data['results'])
                data = requests.get(url=data['next']).json()
            list.extend(data['results'])
            samples = {}
            for item in list:
                samples.setdefault(item['name'], [])
                samples[item["name"]].append(item)

            fuzzy = process.extract(name, samples.keys(), scorer=fuzz.token_set_ratio,
                                    limit=None)
            formatted_results = []
            for name, score in fuzzy:
                for sample in samples[name]:
                    sample.setdefault("relevance", score)
                    formatted_results.append(sample)

            results['results'] = formatted_results

        except:
            traceback.print_exc()
            results["Error"] = "Could not get samples"
    return results
