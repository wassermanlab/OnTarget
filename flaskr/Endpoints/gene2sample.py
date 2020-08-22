import json


def gene2sample(samples, cell_line, treatment, XChrom, YChrom, cancer):
    returnsamples = []
    for sample in samples:
        with open(f'samples/{sample}.json') as f:
            data = json.load(f)
            for item in data:
                for entry in data[item]:
                    if (filterdata(entry['cell_line'], cell_line)) and (filterdata(entry['treatment'], treatment)) and (
                            filterdata(entry['X'], XChrom)) and (filterdata(entry['Y'], YChrom)) and (
                            filterdata(entry['cancer'], cancer)):
                        returnsamples.append(entry)
    return {"results":returnsamples}


def filterdata(data, var):
    if var is not None:
        if data == var:
            return True
        else:
            return False
    else:
        return True
