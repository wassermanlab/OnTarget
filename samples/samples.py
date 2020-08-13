#!/usr/bin/env python

import argparse
import coreapi
import json
import os

usage_msg = """
usage: %s [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
get all samples.

optional arguments:
  -h, --help          show this help message and exit
  -j, --json          output in JSON format
""" % usage_msg

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an
    {argparse} object.
    """

    parser = argparse.ArgumentParser(add_help=False)

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("-j", "--json", action="store_true")
    
    args = parser.parse_args()

    check_args(args)

    return(args)

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

def main():

    # Parse arguments
    args = parse_args()

    # Fuzzy search for samples
    print(get_samples(args.json))

def get_samples(as_json=False):
    """
    e.g. python -m GUD.scripts.get_samples.py"
    """

    try:

        samples = _get_samples()
        samples.sort(key=lambda x: x["name"])

    except:

        raise ValueError("Could not get samples from GUD!!!")

    if as_json:
        return(json.dumps(samples, indent=4))
    else:
        results = set()
        for s in samples:
            results.add(s["name"])
        return("\n".join(r for r in sorted(results)))

def _get_samples():

    # Initialize
    results = []
    client = coreapi.Client()
    codec = coreapi.codecs.CoreJSONCodec()

    # Get the first page
    url = "http://gud.cmmt.ubc.ca:8080/api/v1/hg38/samples"
    response = client.get(url)
    page = json.loads(codec.encode(response))

    # While there are more pages...
    while "next" in page:

        for r in page["results"]:
            results.append(r)

        # Go to the next page
        response = client.get(page["next"])
        page = json.loads(codec.encode(response))

    # Do the last page...
    for r in page["results"]:
        results.append(r)

    return(results)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()
