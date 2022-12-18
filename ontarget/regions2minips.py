#!/usr/bin/env python


from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
import click
from click_option_group import optgroup
import json
import os
import random

from GUD.parsers import ParseUtils

from ontarget import OnTargetUtils
from ontarget.interval2regions import _get_sequence_restriction_sites


class MiniPromoter(object):
    """
    Implements a MiniPromoter object based on the Biopython's Sequence
    Feature object.

    Attributes:
    chrom {str} chromosome of the feature location {FeatureLocation} location
    of the feature on the genome type {str} the specified type of feature (e.g.
    gene, TSS, repeat...) strand {int} on the DNA sequence. \"1\" indicates the
    plus strand; \"-1\" the minus strand; \"0\" for unknown or not applicable
    strand id {str} identifier for the feature qualifiers {dict} qualifiers of
    the feature profile {array} of scores per nucleotide.
    """

    def __init__(
        self,
        chrom,
        starts=[],
        ends=[],
        score=0,
        strand=None,
        feat_type="MiniPromoter",
        feat_ids=["NA"],
        qualifiers=None,
    ):

        self.chrom = chrom
        self._starts = starts
        self._ends = ends
        self.score = score
        self._strand = strand
        self.strand = self.strand_binary
        self.type = feat_type
        self.ids = feat_ids
        self.qualifiers = qualifiers

    @property
    def starts(self):

        return [int(start) for start in self._starts]

    @property
    def starts_1_based(self):

        return [start + 1 for start in self.starts]

    @property
    def ends(self):

        return [int(end) for end in self._ends]

    @property
    def ends_1_based(self):

        return self.ends

    @property
    def strand_binary(self):

        if self._strand == "+":
            return 1
        if self._strand == "-":
            return -1

        return 0

    @property
    def strand_string(self):

        if self.strand == 1:
            return "+"
        if self.strand == -1:
            return "-"

        return "."

    @property
    def gud_strand(self):

        return self._strand

    def __str__(self):

        return "{}\t{},\t{},\t{}\t{}\t{}".\
            format(
                self.chrom,
                ",".join(map(str, self.starts)),  # 0-based for BED format
                ",".join(map(str, self.ends)),
                "+".join(self.ids),
                self.score,
                self.strand_string
            )

    def __repr__(self):

        return "<%s(%s, %s, %s, %s, %s, %s)>" % \
            (
                self.type,
                "chrom={}".format(self.chrom),
                "starts={}".format(",".join(map(str, self.starts))),
                "ends={}".format(",".join(map(str, self.ends))),
                "ids={}".format("+".join(self.ids)),
                "score={}".format(self.score),
                "strand={}".format(self.strand)
            )

    def serialize(self):
        return {
            "chrom": self.chrom,
            "starts": self.starts,
            "ends": self.ends,
            "type": self.type,
            "ids": self.ids,
            "score": self.score,
            "strand": self.strand,
            "qualifiers": self.qualifiers,
        }


CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "json_file",
    type=click.Path(exists=True, resolve_path=True),
)
@optgroup.group("MiniPromoter")
@optgroup.option(
    "--designs",
    help="Max. number of MiniPromoters to design per promoter region.",
    type=int,
    default=OnTargetUtils.max_minip_designs,
    show_default=True,
)
@optgroup.option(
    "--enzyme",
    help="Restriction site to avoid.",
    type=str,
    multiple=True,
)
@optgroup.option(
    "--size",
    help="Max. size of MiniPromoters (in bp).",
    type=int,
    default=OnTargetUtils.max_minip_size,
    show_default=True,
)
@optgroup.option(
    "--tf",
    help="Transcription factor binding site to include.",
    type=str,
    multiple=True,
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

    # Get regulatory regions
    handle = ParseUtils._get_file_handle(args["json_file"])
    regions = json.load(handle)
    handle.close()

    # Get MiniPromoters
    minips = get_minipromoters(regions, args["designs"],
                               set([e.upper() for e in args["enzyme"]]),
                               args["size"],
                               set(tf.upper() for tf in args["tf"]))

    # Write
    json_file = os.path.join(args["output_dir"], "minips.json")
    handle = ParseUtils._get_file_handle(json_file, "w")
    handle.write(json.dumps(minips, indent=4))
    handle.close()


def get_minipromoters(regions, designs=OnTargetUtils.max_minip_designs,
                      enzymes=set(), size=OnTargetUtils.max_minip_size,
                      tfs=set()):
    """
    Function to design MiniPromoters of a max. size from regions identified
    within an interval
    :return: list, GenomicFeatures
    """

    # Initialize
    promoters = []
    enhancers = []
    minips = []

    # Get promoters and enhancers
    for r in regions:
        if enzymes:
            if enzymes.intersection(set(r["qualifiers"]["enzymes"])):
                continue
        if tfs:
            if not tfs.intersection(set(r["qualifiers"]["tfs"])):
                continue
        if r["type"] == "Promoter":
            promoters.append(r)
        else:
            enhancers.append(r)

    # Combine promoters and enhancers
    for promoter in promoters:
        for minip in _get_minipromoters(promoter, enhancers, designs, size):
            minips.append(minip)

    return [minip.serialize() for minip in minips]


def _get_minipromoters(promoter, enhancers,
                       designs=OnTargetUtils.max_minip_designs,
                       size=OnTargetUtils.max_minip_size):

    # Initialize
    minips = []

    # Map enhancer ids to enhancers
    ids2enhancers = {e["id"]:e for e in enhancers}

    # Get promoter strand
    strand = set([tss[1] for tss in promoter["qualifiers"]["TSS"]])
    if len(strand) != 1:
        return minips # i.e., ubiquitous
    strand = strand.pop()

    # Sort enhancers by start position
    enhancers.sort(key=lambda x: x["start"])

    # Get enhancer combinations
    max_size = size - len(promoter["qualifiers"]["sequence"])
    enhancer_combs = [c for c in _enhancer_combinations(enhancers, max_size)]

    # Get MiniPromoters
    random.seed(0)
    for ec in random.sample(enhancer_combs, designs):

        # Get upstream and downstream enhancers
        upstream = []
        downstream = []
        for e in ec:
            if strand == 1: # i.e., forward
                if ids2enhancers[e]["end"] < promoter["start"]:
                    upstream.append(e)
                else:
                    downstream.append(e)
            else: # i.e., reverse
                if ids2enhancers[e]["end"] < promoter["start"]:
                    downstream.append(e)
                else:
                    upstream.append(e)

        # Get MiniPromoter ids
        minip_ids = [promoter["id"]]
        if strand == 1:
            for e in reversed(upstream):
                minip_ids.append(e)
            for e in downstream:
                minip_ids.append(e)
        else:
            for e in upstream:
                minip_ids.append(e)
            for e in reversed(downstream):
                minip_ids.append(e)

        # Get MiniPromoter starts, ends, sequence, score, tfs, etc.
        score = 0
        tfs = set()
        minip_starts = []
        minip_ends = []
        minip_seqs = []
        for minip_id in minip_ids:
            if minip_id in ids2enhancers:
                s = ids2enhancers[minip_id]["qualifiers"]["sequence"]
                if strand == -1:
                    s = str(Seq(s).reverse_complement())
                minip_seqs.insert(0, s)
                score += ids2enhancers[minip_id]["score"] * len(s)
                minip_starts.append(ids2enhancers[minip_id]["start"])
                minip_ends.append(ids2enhancers[minip_id]["end"])
                tfs.update(ids2enhancers[minip_id]["qualifiers"]["tfs"])
            else:
                s = promoter["qualifiers"]["sequence"]
                if strand == -1:
                    s = str(Seq(s).reverse_complement())
                minip_seqs.insert(0, s)
                score += promoter["score"] * len(s)
                minip_starts.append(promoter["start"])
                minip_ends.append(promoter["end"])
                tfs.update(promoter["qualifiers"]["tfs"])
        minip_seq = "".join(minip_seqs)
        tfs = sorted(tfs)
                
        # Get MiniPromoter size
        size = sum(len(s) for s in minip_seqs)

        # Get MiniPromoter score per bp
        score /= size

        # Get MiniPromoter feature
        minip = MiniPromoter(
            chrom=promoter["chrom"],
            starts=minip_starts,
            ends=minip_ends,
            score=score,
            strand="+" if strand == 1 else "-",
            feat_type="MiniPromoter",
            feat_ids=minip_ids,
            qualifiers={
                "enzymes": _get_sequence_restriction_sites(minip_seq),
                "sequence": minip_seq,
                "size": size,
                "tfs": tfs,
            }
        )
        minips.append(minip)

    return minips


def _enhancer_combinations(enhancers, max_size):

    # Initialize
    enhancer_combs = []
    nr_enhancer_combs = []

    # Get enhancer sizes
    sizes = {e["id"]:len(e["qualifiers"]["sequence"]) for e in enhancers}

    # Get enhancer combinations
    for ec in __enhancer_combinations(list(sizes.keys()), sizes, max_size):
        enhancer_combs.append(ec)

    # Sort enhancer combinations by size
    enhancer_combs.sort(key=lambda x: x[-1], reverse=True)

    # Get non-redundant enhancer combinations
    for ec in enhancer_combs:
        is_nr = True
        for nr_ec in nr_enhancer_combs:
            if set(ec[0]).issubset(set(nr_ec)):
                is_nr = False
                break
        if is_nr:
            nr_enhancer_combs.append(ec[0])

    return(nr_enhancer_combs)
    

def __enhancer_combinations(enhancers, sizes, max_size, combination=[],
                            combination_size=0):
    if combination_size > max_size:
        return
    if combination_size > 0:
        yield combination, combination_size
    for i, e in enumerate(enhancers):
        remaining_enhancers = enhancers[i+1:]
        yield from __enhancer_combinations(remaining_enhancers, sizes,
                                           max_size, combination + [e],
                                           combination_size + sizes[e])


if __name__ == "__main__":
    cli()