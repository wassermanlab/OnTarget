#!/usr/bin/env python

import argparse
import coreapi
from fuzzywuzzy import fuzz, process
import json
import os

# Import from OnTarget module
from . import OnTargetUtils

usage_msg = """
usage: %s --name STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
searches one or more samples by string matching.

  --name STR          sample name for string matching

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

    # Mandatory args
    parser.add_argument("--name")

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

    # Check mandatory arguments
    if not args.name:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "argument \"--name\" is required\n"]
        print(": ".join(error))
        exit(0)

def main():

    # Parse arguments
    args = parse_args()

    # Fuzzy search for samples
    print(name_to_sample(args.name, args.json))

def name_to_sample(name, as_json=False):
    """
    e.g. python -m GUD.scripts.name2samples --name "endothelial cells"
    """

    try:

        samples = {}
        for s in _get_samples():
            samples.setdefault(s["name"], [])
            samples[s["name"]].append(s)

    except:

        raise ValueError("Could not get samples from GUD!!!")

    # Fuzzy search
    results = process.extract(name, samples.keys(), scorer=fuzz.token_set_ratio,
        limit=None)

    if as_json:
        formatted_results = []
        for name, score in results:
            for sample in samples[name]:
                sample.setdefault("relevance", score)
                formatted_results.append(sample)
        return(json.dumps(formatted_results, indent=4))
    else:
        return("\n".join("%s\t%s" % (n, s) for n, s in results))

def _get_samples():

    # Initialize
    results = []
    client = coreapi.Client()
    codec = coreapi.codecs.CoreJSONCodec()

    # Get the first page
    response = client.get(os.path.join(OnTargetUtils.gud, "samples"))
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
