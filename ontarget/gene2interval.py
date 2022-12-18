#!/usr/bin/env python


import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import json
import os
import sys

from GUD import GUDUtils
from GUD.ORM import Chrom, Gene
from GUD.ORM.genomic_feature import GenomicFeature
from GUD.parsers import ParseUtils


CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "gene",
    type=str,
)
@click.argument(
    "genome",
    type=click.Choice(["hg19", GUDUtils.db, "mm10"], case_sensitive=True),
)
@optgroup.group(cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "--lim-by-gene",
    help="Limit gene intervals based on the UTRs of the closest " + \
         "upstream and downstream genes.",
    is_flag=True,
)
@optgroup.option(
    "--lim-by-dist",
    help="Limit gene intervals based on distance (i.e., gene body ± `x` kb).",
    type=click.IntRange(0, 1000, clamp=True),
)
@optgroup.group("GUD")
@optgroup.option(
    "--user",
    help="User name",
    type=str,
    default=GUDUtils.user,
    show_default=True,
)
@optgroup.option(
    "--password",
    help="Password.",
    type=str,
    default=GUDUtils.pwd,
    show_default=True,
)
@optgroup.option(
    "--host",
    help="Host name.",
    type=str,
    default=GUDUtils.host,
    show_default=True,
)
@optgroup.option(
    "--port",
    help="Port number.",
    type=int,
    default=GUDUtils.port,
    show_default=True,
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True,
)


def cli(**args):

    # Initialize
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Save exec. parameters as JSON
    json_file = os.path.join(args["output_dir"],
                             f"parameters-{os.path.basename(__file__)}.json")
    handle = ParseUtils._get_file_handle(json_file, "wt")
    handle.write(json.dumps(args, indent=4, sort_keys=True))
    handle.close()

    # Get session
    GUDUtils.user = args["user"]
    GUDUtils.pwd  = args["password"]
    GUDUtils.host = args["host"]
    GUDUtils.port = args["port"]
    GUDUtils.db = args["genome"]
    session = GUDUtils.get_session()

    # Get intervals
    if args["lim_by_gene"]:
        intervals = get_intervals_limit_by_gene(session, args["gene"])
    else:
        intervals = get_intervals_limit_by_distance(session, args["gene"],
                                                    args["lim_by_dist"] * 1000)

    # Write
    json_file = os.path.join(args["output_dir"], "interval.json")
    handle = ParseUtils._get_file_handle(json_file, "w")
    handle.write(f"{intervals}\n")
    handle.close()


def get_intervals_limit_by_gene(session, name):
    """
    Function to get the genomic interval of genes based on the UTRs
    of the closest upstream and downstream genes
    :param session: sqlalchemy Session, session to connect to GUD
    :param name: str, gene name
    :return: list, GenomicFeatures
    """

    # Initialize
    intervals = []

    # Get all chromosomes
    chroms = Chrom.chrom_sizes(session)

    # Get all genes in ncbiRefSeqSelect
    q = Gene.select_by_sources(session, None, ["ncbiRefSeqSelect"])

    # Get gene names
    genes = Gene.select_by_names(session, q, [name])

    # For each gene...
    for g in genes:

        # Get gene feature
        gene = Gene.as_genomic_feature(g)

        # Get chrom, start, end
        chrom = gene.chrom
        start = 0
        end = chroms[chrom]

        # Get all genes in chromosome as gene features
        genes_chrom = [Gene.as_genomic_feature(g2) for g2 in \
            # Gene.select_by_within_location(session, q, chrom, 0, end)]
            Gene.select_by_chrom(session, q, chrom)]

        # Get index of gene
        idx = _binary_search(genes_chrom, 0, len(genes_chrom) - 1, gene)

        # Get start, end
        start = min([gene.start, genes_chrom[idx-1].end]) \
            if idx > 0 else start
        end = max([gene.end, genes_chrom[idx+1].start]) \
            if idx < len(genes_chrom) else end

        # Get interval feature
        feat = GenomicFeature(
            chrom=chrom,
            start=start,
            end=end,
            feat_type="Gene interval",
            feat_id=gene.qualifiers["gene_symbol"],
        )
        intervals.append(feat)

    return json.dumps([i.serialize() for i in intervals], indent=4)


def get_intervals_limit_by_distance(session, name, distance):
    """
    Function to get the genomic interval of genes based on distance
    (i.e., gene body ± `x` kb)
    :param session: sqlalchemy Session, session to connect to GUD
    :param name: str, gene name
    :param distance: int, distance (in kb)  
    :return: list, GenomicFeatures
    """

    # Initialize
    intervals = []

    # Get all genes in ncbiRefSeqSelect
    q = Gene.select_by_sources(session, None, ["ncbiRefSeqSelect"])

    # Get gene names
    genes = Gene.select_by_names(session, q, [name])

    # For each gene...
    for g in genes:

        # Get gene feature
        gene = Gene.as_genomic_feature(g)

        # Get genomic feature
        feat = GenomicFeature(
            chrom=gene.chrom,
            start=gene.start-distance,
            end=gene.end+distance,
            feat_type="Gene interval",
            feat_id="%s " % gene.qualifiers["gene_symbol"] + \
                    f"+/- {int(distance/1000.)} kb"
        )
        intervals.append(feat)

    return json.dumps([i.serialize() for i in intervals], indent=4)


def _binary_search(genes, start, end, gene):

    # Initialize
    middle = (start + end) // 2

    # If gene.start is smaller than middle gene.start, then the feat can only
    # be present in the left sub-array
    if genes[middle].start > gene.start:
        return _binary_search(genes, start, middle - 1, gene)

    # If feat.start is greater than middle gene.start, then the feat can only
    # be present in the right sub-array
    elif genes[middle].start < gene.start:
        return _binary_search(genes, middle + 1, end, gene)

    # Else, this is the feat
    else:
        return middle


if __name__ == "__main__":
    cli()