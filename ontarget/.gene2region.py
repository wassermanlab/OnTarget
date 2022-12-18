#!/usr/bin/env python

import click
from click_option_group import optgroup
import getpass
import json
from multiprocessing import cpu_count
import os
from pybedtools import BedTool

# Import from GUD module
from GUD import GUDUtils

# usage_msg = """
# usage: gene2region.py (--gene [STR ...] | --gene-file FILE)
#                       [-h] [--dummy-dir DIR] [-l] [-o FILE]
#                       [--sample [STR ...] | --sample-file FILE]
#                       [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
# """

# help_msg = """%s
# delimits a region for the given gene(s) based on TAD boundaries,
# UTRs of nearby genes, or genomic distance (in kb).

#   --gene [STR ...]    gene(s) (e.g. "CD19")
#   --gene-file FILE    file containing a list of genes

# optional arguments:
#   -h, --help          show this help message and exit
#   --dummy-dir DIR     dummy directory (default = "/tmp/")
#   -l --limit-by       limit region based on "tad" boundaries
#                       (default), nearby "gene" UTRs, or distance
#   -o FILE             output file (default = stdout)
#   --sample [STR ...]  sample(s) for GUD features (e.g. "B cell")
#   --sample-file FILE  file containing a list of samples for GUD
#                       features

# mysql arguments:
#   -d STR, --db STR    database name (default = "%s")
#   -H STR, --host STR  host name (default = "%s")
#   -p STR, --pwd STR   password (default = ignore this option)
#   -P STR, --port STR  port number (default = %s)
#   -u STR, --user STR  user name (default = "%s")
# """ % \
# (
#     usage_msg,
#     GUDglobals.db_name,
#     GUDglobals.db_host,
#     GUDglobals.db_port,
#     GUDglobals.db_user
# )


CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument("chrom")
@click.argument("start")
@click.argument("end")
@click.option(
    "-c", "--cpu-threads",
    help="Number of CPU threads to use.",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "-d", "--debugging",
    help="Debugging mode.",
    is_flag=True,
)
@click.option(
    "-e", "--evidence",
    help="Genomic evidence.",
    nargs=3,
    type=click.Tuple([click.Path(exists=True, resolve_path=True),
                      click.Choice(["dna_accessibility", "histone_modification", "nascent_transcription", "tf_binding"], case_sensitive=True),
                      click.FloatRange(-1, 1, clamp=True)]),
    multiple=True,
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True,
)

@optgroup.group("GUD")
@optgroup.option(
    "--database",
    help="Database name.",
    type=click.Choice(["hg19", GUDUtils.db, "mm10"], case_sensitive=True),
    default=GUDUtils.db,
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
    "--password",
    help="Password.",
    type=str,
    default=GUDUtils.pwd,
    show_default=True,
)
@optgroup.option(
    "--port",
    help="Port number.",
    type=int,
    default=GUDUtils.port,
    show_default=True,
)
@optgroup.option(
    "--user",
    help="User name",
    type=str,
    default=GUDUtils.user,
    show_default=True,
)

def parse_args():
    """
    This function parses arguments provided via
    the command line and returns an {argparse}
    object.
    """

    parser = argparse.ArgumentParser(
        add_help=False,
    )

    # Mandatory args
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
    optional_group.add_argument(
        "-l", "--limit-by",
        default="tad"
    )
    optional_group.add_argument("-o")
    optional_group.add_argument(
        "--sample", nargs="*"
    )
    optional_group.add_argument("--sample-file")

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
                    "%s\ngene2region.py" % usage_msg,
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
                    "%s\ngene2region.py" % usage_msg,
                    "error",
                    "arguments \"--gene\" \"--gene-file\"",
                    "expected one argument\n"
                ]
            )
        )
        exit(0)

    # Check limit argument
    if not args.limit_by.isdigit():
        if args.limit_by not in ["tad", "gene"]:
            print(": "\
                .join(
                    [
                        "%s\ngene2region.py" % usage_msg,
                        "error",
                        "argument \"--limit-by\"",
                        "invalid choice",
                        "\"%s\" (choose from" % args.limit_by,
                        "\"tad\" \"gene\"; or provide an INT)\n"
                    ]
                )
            )
            exit(0)

    # Check sample arguments
    if args.sample and args.sample_file:
        print(": "\
            .join(
                [
                    "%s\ngene2region.py" % usage_msg,
                    "error",
                    "arguments \"--sample\" \"--sample-file\"",
                    "expected one argument\n"
                ]
            )
        )
        exit(0)

def main(**args):

    pass

    # # Parse arguments
    # args = parse_args()

    # # Initialize
    # dummy_file = os.path.join(
    #     os.path.abspath(args.dummy_dir),
    #     "%s.%s.bed" % (
    #         os.path.basename(__file__),
    #         os.getpid()
    #     )
    # )

    # # Fetch genes
    # genes = []
    # if args.gene:
    #     genes = args.gene
    # else:
    #     gene_file = os.path.abspath(
    #         args.gene_file)
    #     for g in GUDglobals.parse_file(
    #         gene_file):
    #         genes.append(g)

    # # Fetch samples
    # samples = []
    # if args.sample:
    #     samples = args.sample
    # elif args.sample_file:
    #     sample_file = os.path.abspath(
    #         args.sample_file)
    #     for s in GUDglobals.parse_file(
    #         sample_file
    #     ):
    #         samples.append(s)

    # # Establish SQLalchemy session with GUD
    # session = GUDglobals.establish_GUD_session(
    #     args.user,
    #     args.pwd,
    #     args.host,
    #     args.port,
    #     args.db
    # )

    # # For each gene...
    # for gene in genes:
    #     # Get that gene's region
    #     region = get_gene_region(
    #         session,
    #         gene,
    #         samples,
    #         args.limit_by
    #     )
    #     # Write
    #     GUDglobals.write(
    #         dummy_file,
    #         region
    #     )

    # # If output file...
    # if args.o:
    #     shutil.copy(
    #         dummy_file,
    #         os.path.abspath(args.o)
    #     )
    # # ... Else, print on stdout...
    # else:
    #     for line in GUDglobals.parse_file(
    #         dummy_file
    #     ):
    #         GUDglobals.write(None, line)

    # # Delete dummy file
    # os.remove(dummy_file)
    
def get_gene_region(session, gene, samples=[],
    limit_by="tad"):
    """
    Delimits a region for the given gene based
    on TAD boundaries, nearby genes, or +/- N
    kb.
    """

    # Initialize
    chrom = None

    # Get chromosome sizes
    chrom_sizes = Chrom.chrom_sizes(session)

    # Get all genes with the given name
    genes = Gene.select_by_name(session,
        gene, as_genomic_feature=True)
    # If genes are not valid...
    if not genes:
        raise ValueError(
            "Gene \"%s\" is not valid!!!" % gene
        )

    # For each gene...
    for g in genes:

        # Ignore non-standard chroms,
        # scaffolds, etc.
        if g.chrom in chrom_sizes:

            # If first gene...
            if not chrom:
                chrom = g.chrom
                gene_start = chrom_sizes[chrom]
                gene_end = 0
                region_start = 0
                region_end = chrom_sizes[chrom]

            # If gene is on a different chrom...
            if g.chrom != chrom:
                raise ValueError(
                    "Gene \"%s\" is mapped to different chromosomes!!!" % gene
                )

            # If start position is upstream
            # from previous...
            if g.start < gene_start:
                gene_start = g.start

            # If end position is downstream
            # from previous...
            if g.end > gene_end:
                gene_end = g.end

    # If chrom...
    if chrom:

        # If delimit by distance...
        if limit_by.isdigit():
            region_start = gene_start -\
                int(limit_by) * 1000
            region_end = gene_end +\
                int(limit_by) * 1000

        # Instead, if delimit by TADs...
        elif limit_by == "tad":
            region_start, region_end =\
                get_region_coordinates_by_tad(
                    session,
                    gene,
                    chrom,
                    gene_start,
                    gene_end,
                    region_start,
                    region_end,
                    samples
                )

        # Instead, if delimit by genes...
        elif limit_by == "gene":
            region_start, region_end =\
                get_region_coordinates_by_gene(
                    session,
                    gene,
                    chrom,
                    gene_start,
                    gene_end,
                    region_start,
                    region_end
                )

        # Get genomic feature of region
        gene_region = GenomicFeature(
            chrom,
            region_start,
            region_end,
            feat_type = "Region",
            feat_id = "%s" % gene
        )

        return gene_region

    raise ValueError(
        "Gene \"%s\" is not in a valid chromosome!!!" % gene
    )

def get_region_coordinates_by_tad(session,
    gene, chrom, gene_start, gene_end,
    region_start, region_end, samples=[]):
    """
    Delimits a region for the given gene based
    on TAD boundaries.
    """

    # Initialize
    tads = []
    sampleIDs = []

    if samples:
        # Get TADs
        tads = get_tads(
            session,
            chrom,
            gene_start,
            gene_end,
            region_start,
            region_end,
            samples
        )

    if not tads:
        # For each TSS...
        for tss in TSS.select_by_gene(
            session,
            gene,
            as_genomic_feature=True
        ):
            sampleIDs +=\
                tss.qualifiers["sampleIDs"]
        # Get samples
        samples = Sample.select_by_uids(
            session,
            list(set(sampleIDs))
        )
        # Get TADs
        tads = get_tads(
            session,
            chrom,
            gene_start,
            gene_end,
            region_start,
            region_end,
            [s.name for s in samples]
        )

    if not tads:
        # Get TADs
        tads = get_tads(
            session,
            chrom,
            gene_start,
            gene_end,
            region_start,
            region_end
        )

    # For each TAD...
    for t in tads:
        # Get upstream start closest
        # to the gene's start 
        if t.start <= gene_start and \
        t.start > region_start:
            region_start = t.start
        # Get downstream end closest
        # to the gene's end 
        if t.end >= gene_end and \
        t.end < region_end:
            region_end = t.end

    return region_start, region_end

def get_tads(session, chrom, gene_start,
    gene_end, region_start, region_end,
    samples=[]):

    # Initialize
    encompassing_tads = []

    # Get TADs
    tads = TAD.select_by_location(
        session,
        chrom,
        gene_start,
        gene_end,
        samples,
        as_genomic_feature=True
    )

    for t in tads:
        if t.end >= gene_end and \
        t.start <= gene_start:
            encompassing_tads.append(t)

    if encompassing_tads:
        return encompassing_tads

    return tads

def get_region_coordinates_by_gene(session,
    gene, chrom, gene_start, gene_end,
    region_start, region_end):
    """
    Delimits a region for the given gene based
    on nearby genes.
    """

    # Get genes within chromosome
    genes = Gene.select_by_location(
        session,
        chrom,
        region_start,
        region_end,
        as_genomic_feature=True
    )

    # For each gene...
    for g in genes:
        # Skip the gene
        if g.id == gene: continue
        # Get the closest upstream gene end
        # to the gene's start
        if g.end <= gene_start\
        and g.end > region_start:
            region_start = g.end
        # Get the closest downstream gene
        # start to the gene's end
        if g.start >= gene_end\
        and g.start < region_end:
            region_end = g.start

    return region_start, region_end

if __name__ == "__main__":
    main()