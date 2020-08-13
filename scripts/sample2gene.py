#!/usr/bin/env python

import argparse
from functools import partial
from multiprocessing import Pool, cpu_count
import os
import re
import sys

# Import from OnTarget module
from . import OnTargetUtils, human_TFs

usage_msg = """
usage: %s (--sample [STR ...] | --sample-file FILE)
                      [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
searches one or more genes selectively expressed in samples.

  --sample [STR ...]  sample(s) (e.g. "B cell")
  --sample-file FILE  file containing a list of samples

optional arguments:
  -h, --help          show this help message and exit
  -j, --json          output in JSON format
  --threads INT       number of threads to use (default = %s)

expression arguments:
  -a, --all           return TSSs with expression in all samples
                      (default = False)
  -b, --best          return the most selectively expressed TSS
                      for each gene (i.e. the best; default = False)
  --max-tss INT       max. number of TSSs to return (if 0, return
                      all TSSs; default = %s)
  --min-percent FLT   min. percentile of expression for TSS in
                      samples (default = %s)
  --min-tpm FLT       min. expression levels (in TPM) for TSS in
                      input samples (default = %s)
  -t, --tfs           only return TSSs of transcription factors
                      (default = False)
""" % (usage_msg, (cpu_count() - 1), , OnTargetUtils.max_tss,
    OnTargetUtils.min_percent, OnTargetUtils.min_tpm)

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
    parser.add_argument("--sample", nargs="*")
    parser.add_argument("--sample-file")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("-j", "--json", action="store_true")
    optional_group.add_argument("--threads", default=(cpu_count() - 1))

    # Expression args
    exp_group = parser.add_argument_group("expression arguments")
    exp_group.add_argument("-a", "--all", action="store_true")
    exp_group.add_argument("-b", "--best", action="store_true")
    exp_group.add_argument("--max-tss", default=OnTargetUtils.max_tss)    
    exp_group.add_argument("--min-percent", default=OnTargetUtils.min_percent)
    exp_group.add_argument("--min-tpm", default=OnTargetUtils.min_tpm)
    exp_group.add_argument("-t", "--tfs", action="store_true")

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
    if not args.sample and not args.sample_file:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "arguments \"--sample\" \"--sample-file\"", "expected one argument\n"]
        print(": ".join(error))
        exit(0)
    if args.sample and args.sample_file:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "arguments \"--sample\" \"--sample-file\"", "expected one argument\n"]
        print(": ".join(error))
        exit(0)

    # Check "--threads" argument
    try:
        args.threads = int(args.threads)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "argument \"--threads\"", "invalid int value", "\"%s\"\n" % args.threads]
        print(": ".join(error))
        exit(0)

    # Check "--max-tss" argument
    try:
        args.max_tss = int(args.max_tss)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "argument \"--max-tss\"", "invalid int value", "\"%s\"\n" % args.max_tss]
        print(": ".join(error))
        exit(0)

    # Check "--min-percent" argument
    try:
        args.min_percent = float(args.min_percent)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "argument \"--min-percent\"", "invalid float value", "\"%s\"\n" % args.min_percent]
        print(": ".join(error))
        exit(0)

    # Check "--min-tpm" argument
    try:
        args.min_tpm = float(args.min_tpm)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "argument \"--min-tpm\"", "invalid float value", "\"%s\"\n" % args.min_tpm]
        print(": ".join(error))
        exit(0)

def main():

    # Parse arguments
    args = parse_args()

    # Fetch samples
    samples = []
    if args.sample:
        samples = args.sample
    elif args.sample_file:
        with open(args.sample_file, "r") as f:
            for line in f:
                line = line.strip("\n").split("\t")
                samples.append(line[0])
    else:
        raise ValueError("No samples were provided!!!")

    # i.e. endothelial cells
    tmp_answers = [
        "endothelial cell (microvasculature)",
        "endothelial cell (thorax)",
        "endothelial cell (vein)",
        "endothelial cell (artery)",
        "endothelial cell (aorta)",
        "endothelial cell (glomerulus)",
        "endothelial cell (renal glomerulus)",
        "endothelial cell (pulmonary artery)",
        "endothelial cell (lung microvasculature)",
        "endothelial cell (umbilical vein)",
        "endothelial cell (liver sinusoid)",
        "endothelial cell (lymph node)",
        "endothelial cell (brain microvasculature)",
        "endothelial cell (blood vessel dermis)",
        "endothelial cell (microvascular lymphatic vessel dermis)",
        "endothelial progenitor cell (derived from CD14-positive monocyte)"
    ]

    # Identify selectively expressed genes (or TFs)
    g = sample_to_gene(samples, args.json, args.threads, args.all, args.best,
        args.max_tss, args.min_percent, args.min_tpm, args.tfs)

    print(g)

def sample_to_gene(samples, as_json=False, threads=1, all=False, best=False,
    max_tss=OnTargetUtils.max_tss, min_percent=OnTargetUtils.min_percent,
    min_tpm=OnTargetUtils.min_tpm, tfs=False):
    """
    e.g. python -m GUD.scripts.name2samples --name "GM12878"
    """


    print(s)

    # Parse arguments
    global genes
    args = parse_args()

    # Initialize
    dummy_file = "%s.%s.txt" % (os.path.basename(__file__), os.getpid())
    dummy_file = os.path.join(os.path.abspath(args.dummy_dir), dummy_file)



    # Set MySQL options
    GUDUtils.user = args.user
    GUDUtils.pwd = args.pwd
    GUDUtils.host = args.host
    GUDUtils.port = args.port
    GUDUtils.db = args.db

    # Get database name
    db_name = GUDUtils._get_db_name()

    # Get engine/session
    engine, Session = GUDUtils.get_engine_session(db_name)

    # Start a new session
    session = Session()

    # Get selectively expressed TSSs
    TSSs = _get_selectively_expressed_TSSs(session, samples, args.threads, args.all, args.percent, args.tfs, args.tpm)

    # If selectively expressed TSSs...
    if TSSs:

        # Initialize
        genes = set()
        tss_count = 0
        if args.tss <= 0:
            args.tss = len(TSSs)

        # For each TSS...
        for tss in TSSs:

            # Skip if enough TSSs
            if tss_count == args.tss:
                break

            # If best TSS per gene...
            if args.best:

                # Skip if there is a better TSS for this gene
                if tss.qualifiers["gene"] in genes:
                    continue

            # Write
            ParseUtils.write(dummy_file, tss)

            # Add gene and increase count
            genes.add(tss.qualifiers["gene"])
            tss_count += 1

        # If output file...
        if args.o:
            shutil.copy(dummy_file, os.path.abspath(args.o))
        # ... Else, print on stdout...
        else:
            for line in ParseUtils.parse_file(dummy_file):
                ParseUtils.write(None, line)

        # Delete dummy file
        os.remove(dummy_file)

    else:
        raise ValueError("No selectively expressed genes found!!!")

def _get_selectively_expressed_TSSs(session, samples=[], threads=(cpu_count() - 1), exp_in_all_samples=False, percent_exp=0.0, tfs=False, tpm_exp=100):
    """
    Identifies selectively expressed, genic, TSSs.
    """

    # Initialize
    global genes
    global name2uid
    global uid2name
    global sample2tss
    genes = set()
    name2uid = {}
    uid2name = {}
    sample2tss = []
    selectively_expressed_TSSs = []

    # For each sample...
    for feat in Sample.select_by_names(session):
        uid2name.setdefault(feat.uid, feat.name)
        name2uid.setdefault(feat.name, feat.uid)

    # Get genes
    if tfs:
        genes = human_TFs
    else:
        for gene in Gene.get_all_gene_symbols(session).all():
            genes.add(gene[0])

    # For each sample...
    for sample in samples:
        sample2tss.append(set())

    # Get TSS IDs...
    tssIDs = _get_TSS_IDs(session, samples, tpm_exp, exp_in_all_samples)

    # If there are TSS IDs...
    if tssIDs:

        # Select TSSs
        TSSs = TSS.select_by_uids(session, uids=tssIDs).all()

        # Parallelize scan
        pool = Pool(processes=threads)
        parallelized = partial(_filter_TSSs, samples=samples, percent_exp=percent_exp)
        for tss in pool.imap(parallelized, TSSs):
            if tss is not None:
                selectively_expressed_TSSs.append(tss)
        pool.close()
        pool.join()

    # Sort TSSs by percentile expression
    selectively_expressed_TSSs.sort(key=lambda x: float(x.score), reverse=True)

    return(selectively_expressed_TSSs)

def _get_TSS_IDs(session, samples, tpm_exp=100, exp_in_all_samples=False):

    # Initialize
    tssIDs = set()
    query = Expression.select_by_expression(session, tpm_exp)

    # For each TSS...
    for tss in Expression.select_by_samples(session, samples, query).all():

        # Get index
        i = samples.index(tss.Sample.name)
        sample2tss[i].add(tss.Expression.tssID)
        tssIDs.add(tss.Expression.tssID)

        # Enables search for cell line-specific TSSs
        uid2name.setdefault(tss.Sample.uid, tss.Sample.name)
        name2uid.setdefault(tss.Sample.name, tss.Sample.uid)

    # If required expression in all samples...
    if exp_in_all_samples:
        return(list(set.intersection(*sample2tss)))
    # ... Else...
    else:
        return(list(tssIDs))

def _filter_TSSs(tss, samples, percent_exp=0.0):
    """
    Returns selectively expressed TSSs.
    """

    # Initialize
    tss = TSS.as_genomic_feature(tss)

    # If gene TSS...
    if tss.qualifiers["gene"] in genes:

        # Initialize
        bexp = 0.0 # background expression levels
        fexp = 0.0 # foreground expression levels
        tss.qualifiers.setdefault("expression", {})

        # Compile regexps
        int_regexp = re.compile("\d+")
        float_regexp = re.compile("([0-9]*[.])?[0-9]+")

        # Get sample IDs, avg. expression levels
        sampleIDs = list(map(float, re.findall(int_regexp, str(tss.qualifiers["sampleIDs"]))))
        avg_expression_levels = list(map(float, re.findall(float_regexp, str(tss.qualifiers["avg_expression_levels"]))))

        # For each sampleID
        for i in range(len(sampleIDs)):

            # Initialize
            sampleID = sampleIDs[i]
            avg_expression = avg_expression_levels[i]

            # If sample...
            if sampleID in uid2name:

                # Foreground expression
                if uid2name[sampleID] in samples:
                    fexp += avg_expression

                # Background expression
                bexp += avg_expression

                # Sample expression
                tss.qualifiers["expression"].setdefault(uid2name[sampleID], [avg_expression, None])

        # If expressed...
        if bexp > 0:

            # If selectively expressed...
            if (fexp * 100.0) / bexp >= percent_exp:

                # Initialize 
                tss.id = "p{}@{}".format(tss.qualifiers["tss"], tss.qualifiers["gene"])
                tss.score = "%.3f" % round((fexp * 100.0) / bexp, 3)

                # For each sample...
                for sample in tss.qualifiers["expression"]:
                    tss.qualifiers["expression"][sample][1] = round(tss.qualifiers["expression"][sample][0] * 100.0 / bexp, 3)

                return(tss)

    return(None)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()