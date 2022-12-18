#!/usr/bin/env python

import argparse
import os
import re
import shutil
import sys

# Import from GUD
from GUD import GUDglobals
from GUD.ORM.expression import Expression
from GUD.ORM.gene import Gene
from GUD.ORM.sample import Sample
from GUD.ORM.tss import TSS

usage_msg = """
usage: gene2sample.py (--gene STR ... | --gene-file FILE)
                      [-h] [--dummy-dir DIR] [-o FILE]
                      [--percent FLT] [--tpm FLT]
                      [--sample [STR ...] | --sample-file FILE]
                      [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
"""

help_msg = """%s
identifies one or more samples in which genes are expressed.

  --gene [STR ...]    gene(s) (e.g. "CD19")
  --gene-file FILE    file containing a list of genes

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -o FILE             output file (default = stdout)

expression arguments:
  --percent FLT       min. percentile of expression for TSS in
                      output samples (default = %s)
  --sample [STR ...]  sample(s) (e.g. "B cell")
  --sample-file FILE  file containing a list of samples
  --tpm FLT           min. expression levels (in TPM) for TSS in
                      output samples (default = %s)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "%s")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
  -u STR, --user STR  user name (default = "%s")
""" % \
(
    usage_msg,
    GUDglobals.min_percent,
    GUDglobals.min_tpm,
    GUDglobals.db_name,
    GUDglobals.db_host,
    GUDglobals.db_port,
    GUDglobals.db_user
)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via
    the command line and returns an {argparse}
    object.
    """

    parser = argparse.ArgumentParser(
        add_help=False,
    )

    # Mandatory arguments
    parser.add_argument(
        "--gene", nargs="*"
    )
    parser.add_argument("--gene-file")

    # Optional args
    optional_group = parser.add_argument_group(
        "optional arguments"
    )
    optional_group.add_argument(
        "-h", "--help",
        action="store_true"
    )
    optional_group.add_argument(
        "--dummy-dir",
        default="/tmp/"
    )
    optional_group.add_argument("-o")

    # Expression args
    exp_group = parser.add_argument_group(
        "expression arguments"
    )
    exp_group.add_argument(
        "--percent",
        default=GUDglobals.min_percent
    )
    exp_group.add_argument(
        "--sample", nargs="*"
    )
    exp_group.add_argument("--sample-file")
    exp_group.add_argument(
        "--tpm",
        default=GUDglobals.min_tpm
    )

    # MySQL args
    mysql_group = parser.add_argument_group(
        "mysql arguments"
    )
    mysql_group.add_argument(
        "-d", "--db",
        default=GUDglobals.db_name,
    )
    mysql_group.add_argument(
        "-H", "--host",
        default=GUDglobals.db_host
    )
    mysql_group.add_argument("-p", "--pwd")
    mysql_group.add_argument(
        "-P", "--port",
        default=GUDglobals.db_port
    )
    mysql_group.add_argument(
        "-u", "--user",
        default=GUDglobals.db_user
    )

    args = parser.parse_args()

    check_args(args)

    return args

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Print help
    if args.help:
        print(help_msg)
        exit(0)
    
    # Check mandatory arguments
    if (
        not args.gene and \
        not args.gene_file
    ):
        print(": "\
            .join(
                [
                    "%s\ngene2sample.py" % usage_msg,
                    "error",
                    "one of the arguments \"--gene\" \"--gene-file\" is required\n"
                ]
            )
        )
        exit(0)

    if args.gene and args.gene_file:
        print(": "\
            .join(
                [
                    "%s\ngene2sample.py" % usage_msg,
                    "error",
                    "arguments \"--gene\" \"--gene-file\"",
                    "expected one argument\n"
                ]
            )
        )
        exit(0)

    # Check for invalid percent
    try:
        args.percent = float(args.percent)
    except:
        print(": "\
            .join(
                [
                    "%s\nsample2gene.py" % usage_msg,
                    "error",
                    "argument \"--percent\"",
                    "invalid float value",
                    "\"%s\"\n" % args.percent
                ]
            )
        )
        exit(0)

    # Check for invalid TPM
    try:
        args.tpm = float(args.tpm)
    except:
        print(": "\
            .join(
                [
                    "%s\ngene2sample.py" % usage_msg,
                    "error",
                    "argument \"--tpm\"",
                    "invalid float value",
                    "\"%s\"\n" % args.tpm
                ]
            )
        )
        exit(0)

    if args.sample and args.sample_file:
        print(": "\
            .join(
                [
                    "%s\nsample2gene.py" % usage_msg,
                    "error",
                    "arguments \"--sample\" \"--sample-file\"",
                    "expected one argument\n"
                ]
            )
        )
        exit(0)

def main():

    # Parse arguments
    args = parse_args()
    

    # Initialize
    dummy_file = os.path.join(
        os.path.abspath(args.dummy_dir),
        "%s.%s.txt" % (os.path.basename(__file__),
        os.getpid()))

    # Fetch genes
    genes = []
    if args.gene:
        genes = args.gene
    elif args.gene_file:
        gene_file = os.path.abspath(args.gene_file)
        for gene in GUDglobals.parse_file(gene_file):
            genes.append(gene)
    else:
        raise ValueError(
            "No genes were provided!!!"
        )

    # Fetch samples
    allowed_samples = set()
    if args.sample:
        allowed_samples = set(args.sample)
    elif args.sample_file:
        sample_file = os.path.abspath(args.sample_file)
        for sample in GUDglobals.parse_file(sample_file):
            allowed_samples.add(sample)

    # Establish SQLalchemy session with GUD
    session = GUDglobals.establish_GUD_session(
        args.user,
        args.pwd,
        args.host,
        args.port,
        args.db
    )

    # Get samples from GUD
    samples = {}
    for s in Sample.select_by_names(session):
        samples.setdefault(s.uid, s.name)
    if not allowed_samples:
        allowed_samples = set(samples.values())

    # For each gene...
    for gene in genes:

        # For each TSS...
        for tss in TSS.select_by_gene(
            session,
            gene,
            as_genomic_feature=True
        ):

            # Write
            GUDglobals.write(
                dummy_file,
                ">%s;TSS%s" % (gene, tss.qualifiers["tss"])
            )

            # Get TSS expression
            expression = get_tss_expression(
                session,
                tss
            )

            # For each sample...
            for s in sorted(
                expression,
                key=lambda x: expression[x],
                reverse=True
            ):

                # Get percent expression
                percent_exp = expression[s] * 100
                percent_exp /= sum([
                    expression[s] for s in expression \
                        if s in samples
                ])
                # If expressed...
                if (
                    percent_exp >= args.percent and \
                    expression[s] >= args.tpm and \
                    s in samples
                ):
                    # Skip if not allowed sample
                    if samples[s] not in allowed_samples:
                        continue
                    # Write
                    GUDglobals.write(
                        dummy_file,
                        "%s\t%s\t%s" % (
                            samples[s],
                            expression[s],
                            "%.3f%%" % percent_exp
                        )
                    )

            # Write
            GUDglobals.write(
                dummy_file,
                "//"
            )

    # If output file...
    if args.o:
        shutil.copy(
            dummy_file,
            os.path.abspath(args.o)
        )
    # ... Else, print on stdout...
    else:
        for line in GUDglobals.parse_file(dummy_file):
            GUDglobals.write(None, line)

    # Delete dummy file
    os.remove(dummy_file)

def get_tss_expression(session, tss):
    """
    Identifies samples where TSS is expressed.
    """

    # Initialize
    expression = {}
    isfloat = re.compile("\d+(\.\d+)?")

    # Get sample IDs, avg. expression levels
    # and indices
    sampleIDs = str(
        tss.qualifiers["sampleIDs"]
    ).split(",")
    avg_exps = str(
        tss.qualifiers["avg_expression_levels"]
    ).split(",")
    idxs = list(reversed(range(len(sampleIDs))))

    # For each sample ID...
    for i in idxs:
        if (
            sampleIDs[i].isdigit() and \
            isfloat.match(avg_exps[i])
        ):
            sampleIDs[i] = int(sampleIDs[i])
            avg_exps[i] = float(avg_exps[i])
        else:
            sampleIDs.pop(i)
            avg_exps.pop(i)

    # For each sampleID
    for i in range(len(sampleIDs)):
        # Add sample expression
        expression.setdefault(
            sampleIDs[i],
            avg_exps[i]
        )

    return expression

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()