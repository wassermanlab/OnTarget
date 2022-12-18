#!/usr/bin/env python

#!/usr/bin/env python

import argparse
import coreapi
from fuzzywuzzy import fuzz, process
import json
import os

# Import from OnTarget module
from . import OnTargetUtils

usage_msg = """
usage: %s (--gene [STR ...] | --gene-file FILE)
                      [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
delimits a region for the given gene(s) based on TAD boundaries,
UTRs of nearby genes, or genomic distance (in kb).

  --gene [STR ...]    gene(s) (e.g. "CD19")
  --gene-file FILE    file containing a list of genes

optional arguments:
  -h, --help          show this help message and exit
  -j, --json          output in JSON format
  -l, --limit         limit region based on "tad" boundaries
                      (default), nearby "gene" UTRs, or distance
  --sample [STR ...]  sample(s) for GUD features (e.g. "B cell")
  --sample-file FILE  file containing a list of samples for GUD
                      features
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
    parser.add_argument("--gene", nargs="*")
    parser.add_argument("--gene-file")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("-j", "--json", action="store_true")
    optional_group.add_argument("-l", "--limit", default="tad")
    optional_group.add_argument("--sample", nargs="*")
    optional_group.add_argument("--sample-file")

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
    if (not args.gene and not args.gene_file) or (args.gene and args.gene_file):
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "arguments \"--gene\" \"--gene-file\"", "expected one argument\n"]
        print(": ".join(error))
        exit(0)

    # Check limit argument
    if not args.limit.isdigit():
        if args.limit not in ["tad", "gene"]:
            error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
                "arguments \"--limit\"", "invalid choice", "\"%s\" (choose from" % args.limit,
                "\"tad\" \"gene\"; or provide an INT)\n"]
            print(": ".join(error))
            exit(0)

    # Check sample arguments
    if args.sample and args.sample_file:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "arguments \"--sample\" \"--sample-file\"", "expected one argument\n"]
        print(": ".join(error))
        exit(0)

def main():

    # Parse arguments
    args = parse_args()

    # Fetch genes
    genes = []
    if args.gene:
        genes = args.gene
    elif args.gene_file:
        with open(args.gene_file, "r") as f:
            for g in f:
                g = g.strip("\n")
                genes.append(g)

    # Fetch samples
    samples = []
    if args.sample:
        samples = args.sample
    elif args.sample_file:
        with open(args.sample_file, "r") as f:
            for s in f:
                s = s.strip("\n")
                samples.append(s)

    for g in genes:

        # Delimit region for input gene(s)
        print(gene_to_region(g, args.json, args.limit, samples))
  
def gene_to_region(gene, as_json=False, limit="tad", samples=[]):
    """
    Delimits a region for the given gene based on TAD boundaries,
    UTRs of nearby genes, or genomic distance (in kb).
    """

    try:

        chrom_sizes = {}
        for c in _get_chroms():
            chrom_sizes.setdefault(c["chrom"], c["size"])

    except:

        raise ValueError("Could not get chromosomes from GUD!!!")

    print(chrom_sizes)
    exit(0)

#     # Get all genes with the given name
#     genes = Gene.select_by_name(session,
#         gene, as_genomic_feature=True)
#     # If genes are not valid...
#     if not genes:
#         raise ValueError(
#             "Gene \"%s\" is not valid!!!" % gene
#         )

#     # For each gene...
#     for g in genes:

#         # Ignore non-standard chroms,
#         # scaffolds, etc.
#         if g.chrom in chrom_sizes:

#             # If first gene...
#             if not chrom:
#                 chrom = g.chrom
#                 gene_start = chrom_sizes[chrom]
#                 gene_end = 0
#                 region_start = 0
#                 region_end = chrom_sizes[chrom]

#             # If gene is on a different chrom...
#             if g.chrom != chrom:
#                 raise ValueError(
#                     "Gene \"%s\" is mapped to different chromosomes!!!" % gene
#                 )

#             # If start position is upstream
#             # from previous...
#             if g.start < gene_start:
#                 gene_start = g.start

#             # If end position is downstream
#             # from previous...
#             if g.end > gene_end:
#                 gene_end = g.end

#     # If chrom...
#     if chrom:

#         # If delimit by distance...
#         if limit_by.isdigit():
#             region_start = gene_start -\
#                 int(limit_by) * 1000
#             region_end = gene_end +\
#                 int(limit_by) * 1000

#         # Instead, if delimit by TADs...
#         elif limit_by == "tad":
#             region_start, region_end =\
#                 get_region_coordinates_by_tad(
#                     session,
#                     gene,
#                     chrom,
#                     gene_start,
#                     gene_end,
#                     region_start,
#                     region_end,
#                     samples
#                 )

#         # Instead, if delimit by genes...
#         elif limit_by == "gene":
#             region_start, region_end =\
#                 get_region_coordinates_by_gene(
#                     session,
#                     gene,
#                     chrom,
#                     gene_start,
#                     gene_end,
#                     region_start,
#                     region_end
#                 )

#         # Get genomic feature of region
#         gene_region = GenomicFeature(
#             chrom,
#             region_start,
#             region_end,
#             feat_type = "Region",
#             feat_id = "%s" % gene
#         )

#         return gene_region

#     raise ValueError(
#         "Gene \"%s\" is not in a valid chromosome!!!" % gene
#     )

# def get_region_coordinates_by_tad(session,
#     gene, chrom, gene_start, gene_end,
#     region_start, region_end, samples=[]):
#     """
#     Delimits a region for the given gene based
#     on TAD boundaries.
#     """

#     # Initialize
#     tads = []
#     sampleIDs = []

#     if samples:
#         # Get TADs
#         tads = get_tads(
#             session,
#             chrom,
#             gene_start,
#             gene_end,
#             region_start,
#             region_end,
#             samples
#         )

#     if not tads:
#         # For each TSS...
#         for tss in TSS.select_by_gene(
#             session,
#             gene,
#             as_genomic_feature=True
#         ):
#             sampleIDs +=\
#                 tss.qualifiers["sampleIDs"]
#         # Get samples
#         samples = Sample.select_by_uids(
#             session,
#             list(set(sampleIDs))
#         )
#         # Get TADs
#         tads = get_tads(
#             session,
#             chrom,
#             gene_start,
#             gene_end,
#             region_start,
#             region_end,
#             [s.name for s in samples]
#         )

#     if not tads:
#         # Get TADs
#         tads = get_tads(
#             session,
#             chrom,
#             gene_start,
#             gene_end,
#             region_start,
#             region_end
#         )

#     # For each TAD...
#     for t in tads:
#         # Get upstream start closest
#         # to the gene's start 
#         if t.start <= gene_start and \
#         t.start > region_start:
#             region_start = t.start
#         # Get downstream end closest
#         # to the gene's end 
#         if t.end >= gene_end and \
#         t.end < region_end:
#             region_end = t.end

#     return region_start, region_end

# def get_tads(session, chrom, gene_start,
#     gene_end, region_start, region_end,
#     samples=[]):

#     # Initialize
#     encompassing_tads = []

#     # Get TADs
#     tads = TAD.select_by_location(
#         session,
#         chrom,
#         gene_start,
#         gene_end,
#         samples,
#         as_genomic_feature=True
#     )

#     for t in tads:
#         if t.end >= gene_end and \
#         t.start <= gene_start:
#             encompassing_tads.append(t)

#     if encompassing_tads:
#         return encompassing_tads

#     return tads

# def get_region_coordinates_by_gene(session,
#     gene, chrom, gene_start, gene_end,
#     region_start, region_end):
#     """
#     Delimits a region for the given gene based
#     on nearby genes.
#     """

#     # Get genes within chromosome
#     genes = Gene.select_by_location(
#         session,
#         chrom,
#         region_start,
#         region_end,
#         as_genomic_feature=True
#     )

#     # For each gene...
#     for g in genes:
#         # Skip the gene
#         if g.id == gene: continue
#         # Get the closest upstream gene end
#         # to the gene's start
#         if g.end <= gene_start\
#         and g.end > region_start:
#             region_start = g.end
#         # Get the closest downstream gene
#         # start to the gene's end
#         if g.start >= gene_end\
#         and g.start < region_end:
#             region_end = g.start

#     return region_start, region_end

def _get_chroms():

    # Initialize
    results = []
    client = coreapi.Client()
    codec = coreapi.codecs.CoreJSONCodec()

    # Get the first page
    response = client.get(os.path.join(OnTargetUtils.gud, "chroms"))
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
