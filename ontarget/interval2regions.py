#!/usr/bin/env python


from Bio.Restriction import *
from Bio.Seq import Seq
import click
from click_option_group import optgroup
import json
import numpy as np
import os
import pandas as pd
from pybedtools import BedTool
import re
import requests
from sklearn.preprocessing import MinMaxScaler
import subprocess as sp
import sys

from GUD import GUDUtils
from GUD.ORM import Conservation, Gene, Region, RepeatMask
from GUD.ORM.genomic_feature import GenomicFeature
from GUD.parsers import ParseUtils

from ontarget import OnTargetUtils


CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "chrom",
    type=str,
)
@click.argument(
    "start",
    type=int,
)
@click.argument(
    "end",
    type=int,
)
@click.argument(
    "genome",
    type=click.Choice(["hg19", GUDUtils.db, "mm10"], case_sensitive=True),
)
@optgroup.group("Regulatory regions")
@optgroup.option(
    "--evidence",
    help="Genomic evidence (i.e., `file weight`).",
    nargs=2,
    type=click.Tuple([click.Path(exists=True, resolve_path=True),
                      click.FloatRange(-1, 1, clamp=True)]),
    multiple=True,
)
@optgroup.option(
    "--liftover",
    help="LiftOver regulatory regions to genome.",
    type=click.Choice(["hg19", GUDUtils.db, "mm10"], case_sensitive=True),
)
@optgroup.option(
    "--region-length",
    help="Min. length for regulatory regions.",
    type=int,
    default=OnTargetUtils.get_min_length(),
    show_default=True,
)
@optgroup.option(
    "--region-score",
    help="Min. regulatory region score. If `x` > 1, use score of " +
         "top `x` percentile bp.",
    type=click.FloatRange(0., 100., clamp=True),
    default=OnTargetUtils.get_min_score(),
    show_default=True,
)
@optgroup.group("Conservation")
@optgroup.option(
    "--conservation",
    help="Use conservation from multiz alignments.",
    is_flag=True,
)
@optgroup.option(
    "--cons-score",
    help="Min. conserved region score.",
    type=click.FloatRange(0., 1., clamp=True),
    default=OnTargetUtils.get_min_score("conserved region"),
    show_default=True,
)
@optgroup.option(
    "--cons-length",
    help="Min. length for conserved regions.",
    type=int,
    default=OnTargetUtils.get_min_length("conserved region"),
    show_default=True,
)
@optgroup.group("Masking")
@optgroup.option(
    "--mask-exons",
    help="Mask coding exons.",
    is_flag=True,
)
@optgroup.option(
    "--mask-repeats",
    help="Mask repeating elements (e.g., SINEs, LINEs, etc.) from rmsk.",
    is_flag=True,
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
    "--passwd",
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
    "-d", "--dummy-dir",
    help="Dummy directory.",
    type=click.Path(resolve_path=True),
    default=OnTargetUtils.get_dir("dummy"),
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
    session = OnTargetUtils.get_gud_session(args["genome"], args["user"],
                                            args["passwd"], args["host"],
                                            args["port"])

    # Get regulatory regions
    regions = get_regions(session, args["chrom"], args["start"],
                          args["end"], args["genome"], args["evidence"],
                          args["liftover"], args["region_length"],
                          args["region_score"], args["conservation"],
                          args["cons_score"], args["cons_length"],
                          args["mask_exons"], args["mask_repeats"],
                          args["dummy_dir"])

    # Write
    bed_file = os.path.join(args["output_dir"], "regions.bed")
    OnTargetUtils.write_bed(regions, bed_file)
    # fasta_file = os.path.join(args["output_dir"], "regions.fa")
    # OnTargetUtils.write_fasta(regions, fasta_file)
    json_file = os.path.join(args["output_dir"], "regions.json")
    OnTargetUtils.write_json(regions, json_file)


def get_regions(session, chrom, start, end, genome, evidence=[],
                liftover=None, region_length=OnTargetUtils.get_min_length(),
                region_score=OnTargetUtils.get_min_score(),
                use_conservation=False,
                cons_score=OnTargetUtils.get_min_score("conserved region"),
                cons_length=OnTargetUtils.get_min_length("conserved region"),
                mask_exons=False, mask_repeats=False, dummy_dir="/tmp/"):
    """
    Function to get the regulatory regions within an interval based on 
    genomic evidence
    :param session: SQLAlchemy Session, session to connect to GUD
    :param chrom: str, chromosome
    :param start: int, start coordinate (0-based)
    :param end: int, end coordinate
    :param genome: str, genome assembly
    :param evidence: list, genomic evidence files and weights as sublists
    :param liftover: str, genome assembly to liftOver to
    :param region_length: int, min. regulatory region length
    :param region_score: float, min. regulatory region score 
    :param use_conservation: bool, whether or not to use conservation
    :param cons_score: int, min. conserved region length
    :param cons_length: float, min. conserved region score
    :param mask_exons: bool, whether or not to mask exons
    :param mask_repeats: bool, whether or not to mask repeats
    :param dummy_dir: str, path to dummy directory
    :return: list, regulatory regions as serialied GenomicFeature
    """

    # Get all transcription start sites in interval
    tsss = _get_interval_transcription_start_sites(session, chrom, start, end)

    # Get all coding exons in interval
    exons = _get_interval_coding_exons(session, chrom, start, end)

    # Get evidence profiles
    evidence_feats = _get_evidence_features(chrom, start, end, evidence)
    profiles = [features_to_profile(start, end, f) for f in evidence_feats]

    # Get conservation profiles
    if use_conservation:
        source_name = "multiz60way" if genome == "mm10" else "multiz100way"
        cons_feats = _get_conserved_features(session, chrom, start, end,
                                             genome)
        cons_profile = features_to_profile(start, end, cons_feats)
        profiles.append(cons_profile)
        cons_regions = profile_to_regions(chrom, start, end, genome,
                                          cons_profile, cons_score, cons_length,
                                          exons=exons, stitch=True,
                                          source=source_name)
        cons_regions_profile = features_to_profile(start, end, cons_regions)
        cons_regions_profile[cons_regions_profile > 0.] = 1. # set to 1.
        profiles.append(cons_regions_profile)

    # Build a profile matrix
    profile_matrix = np.array(profiles)

    # Sum each column in the profile matrix to get a vector of the sum
    # of the individual scores at each position
    reg_profile = np.sum(profile_matrix, axis=0)

    # Mask exons and repeats
    mask_profile = np.ones(end - start)
    if mask_exons:
        exon_profile = features_to_profile(start, end, exons)
        mask_profile[exon_profile == 1.] = 0.
    if mask_repeats:
        rmsk_feats = _get_rmsk_features(session, chrom, start, end, genome)
        rmsk_profile = features_to_profile(start, end, rmsk_feats)
        mask_profile[rmsk_profile == 1.] = 0.
    reg_profile *= mask_profile

    # Normalize
    scaler = MinMaxScaler()
    reg_profile_norm = scaler.fit_transform(reg_profile.reshape(-1, 1))\
                             .flatten()

    # Use the top score percentile to determine the scoring threshold used to
    # define regulatory regions
    if region_score >= 1:
        sorted_profile = np.sort(reg_profile_norm)
        top_pct_pos = int(sorted_profile.size - sorted_profile.size * \
                          region_score / 100)
        region_score = max([sorted_profile[top_pct_pos],
                            min(sorted_profile[sorted_profile > 0])])
    
    # Compute regulatory regions
    reg_regions = profile_to_regions(chrom, start, end, genome,
                                     reg_profile_norm, region_score,
                                     region_length, exons=exons,
                                     feat_type="Enhancer",
                                     feat_id_prefix="RR",
                                     source="OnTarget", reset_count=True)

    # Get promoters
    for i in range(len(reg_regions)):
        for tss in tsss:
            if reg_regions[i].overlaps(tss):
                reg_regions[i].type = "Promoter"
                reg_regions[i].qualifiers.setdefault("TSS", set())
                reg_regions[i].qualifiers["TSS"].add((tss.id, tss.strand))
        if reg_regions[i].type == "Promoter":
            reg_regions[i].qualifiers["TSS"] = [tss for tss in \
                reg_regions[i].qualifiers["TSS"]]

    # LiftOver
    if liftover and genome != liftover:
        reg_regions = _liftover_regions(reg_regions, genome, liftover,
                                        dummy_dir)
        genome = liftover

    # Get sequences and restriction sites
    for i in range(len(reg_regions)):
        s = str(_get_interval_sequence(reg_regions[i].chrom,
                                       reg_regions[i].start,
                                       reg_regions[i].end, genome)).upper()
        rs = _get_sequence_restriction_sites(s)
        reg_regions[i].qualifiers.setdefault("sequence", s)
        reg_regions[i].qualifiers.setdefault("enzymes", rs)

    # Get TFs
    for i in range(len(reg_regions)):
        tfs = _get_interval_tfs(reg_regions[i].chrom, reg_regions[i].start,
                                reg_regions[i].end, genome)
        reg_regions[i].qualifiers.setdefault("tfs", list(tfs))

    return [rr.serialize() for rr in reg_regions]


def _get_interval_transcription_start_sites(session, chrom, start, end):

    # Initialize
    tsss = []
    seen_tsss = set()

    # Get all genes
    q = Gene.make_query(session, None)

    # Get genes in interval
    genes = Gene.select_by_location(session, q, chrom, start, end,
                                    location="overlapping")

    # For each gene...
    for g in genes:

        # Get gene feature
        gene = Gene.as_genomic_feature(g)

        # TSS start, end
        if gene.strand_string == "+":
            start = gene.start
            end = start + 1
        else:
            end = gene.end
            start = end - 1

        # Skip
        if (start, end, gene.qualifiers["gene_symbol"]) in seen_tsss:
            continue
        else:
            seen_tsss.add((start, end, gene.qualifiers["gene_symbol"]))

        # Get genomic feature
        feat = GenomicFeature(
            chrom=gene.chrom,
            start=start,
            end=end,
            score=1.,
            strand=gene.strand_string,
            feat_type="TSS",
            feat_id=gene.qualifiers["gene_symbol"],
        )
        tsss.append(feat)
    
    return tsss


def _get_interval_coding_exons(session, chrom, start, end):

    # Initialize
    exons = []
    seen_exons = set()

    # Get all genes
    q = Gene.make_query(session, None)

    # Get genes in interval
    genes = Gene.select_by_location(session, q, chrom, start, end,
                                    location="overlapping")

    # For each gene...
    for g in genes:

        # Get gene feature
        gene = Gene.as_genomic_feature(g)

        # For each exon...
        for i, (start, end) in enumerate(zip(gene.qualifiers["exon_starts"],
                                             gene.qualifiers["exon_ends"])):

            # Skip
            if (start, end, gene.qualifiers["gene_symbol"]) in seen_exons:
                continue
            else:
                seen_exons.add((start, end, gene.qualifiers["gene_symbol"]))
            if int(start) >= gene.qualifiers["coding_end"]:
                continue
            if int(end) <= gene.qualifiers["coding_start"]:
                continue

            # Fix coding start, end
            if start <= gene.qualifiers["coding_start"] and \
               end >= gene.qualifiers["coding_start"]:
               start = gene.qualifiers["coding_start"]
            if end >= gene.qualifiers["coding_end"] and \
               start <= gene.qualifiers["coding_end"]:
               end = gene.qualifiers["coding_end"]
    
            # Get genomic feature
            feat = GenomicFeature(
                chrom=gene.chrom,
                start=int(start),
                end=int(end),
                score=1.,
                strand=gene.strand_string,
                feat_type="Exon",
                feat_id=gene.qualifiers["gene_symbol"],
            )
            exons.append(feat)

    return exons


def _get_interval_sequence(chrom, start, end, genome):
    """
    Function to get the sequence of an interval from UCSC
    :param chrom: str, chromosome name
    :param start: int, start coordinate (0-based)
    :param end: int, end coordinate
    :param genome: str, genome assembly in UCSC Genome Browser
    :return: Bio Seq, sequence
    """

    # Get URL
    url = "https://api.genome.ucsc.edu/getData/" + \
         f"sequence?genome={genome};chrom=chr{chrom};start={start};end={end}"

    # Get response
    response = requests.get(url)

    return Seq(response.json()["dna"])


def _get_sequence_restriction_sites(sequence):
    return [str(re).upper() for re in \
            Analysis(AllEnzymes, Seq(sequence)).with_sites()]


def _get_conserved_features(session, chrom, start, end, genome):
    """
    Function to get the sequence of an interval
    :param session: SQLAlchemy Session, session to connect to GUD
    :return: Bio Seq, sequence
    """

    # Initialize
    cons_feats = []

    # Get all conserved elements in multiz alignments
    source_name = "multiz60way" if genome == "mm10" else "multiz100way"
    q = Conservation.select_by_sources(session, None, [source_name])

    # Get conserved elements in interval
    cons_elems = Conservation.select_by_location(session, q, chrom, start, end,
                                                 location="overlapping")

    # For each conserved element...
    for ce in cons_elems:

        # Get conserved feature
        cons_feat = Conservation.as_genomic_feature(ce)
        cons_feats.append(cons_feat)

    return cons_feats


def _get_rmsk_features(session, chrom, start, end, genome):
    """
    Function to get the sequence of an interval
    :param session: SQLAlchemy Session, session to connect to GUD
    :return: Bio Seq, sequence
    """

    # Initialize
    rmsk_feats = []

    # Get rmsk elements in interval
    q = RepeatMask.make_query(session, None)
    rmsk_elems = RepeatMask.select_by_location(session, q, chrom, start, end,
                                               location="overlapping")

    # For each rmsk element...
    for re in rmsk_elems:

        # Get rmsk feature
        rmsk_feat = RepeatMask.as_genomic_feature(re)
        rmsk_feat.score = 1.
        rmsk_feats.append(rmsk_feat)

    return rmsk_feats


def _get_evidence_features(chrom, start, end, evidence=[]):

    # Initialize
    # evidence_feats = {}
    evidence_feats = []

    # For each evidence...
    # for bed_file, evidence, weight in evidence:
    for bed_file, weight in evidence:


        # Initialize
        feat_type = os.path.basename(bed_file)
        # evidence_feats.setdefault(evidence, [])
        # evidence_feats[evidence].append([])
        evidence_feats.append([])

        # Open file as DataFrame
        df = pd.read_table(bed_file, header=None, usecols=[0, 1, 2])

        # Select overlapping intervals
        df = df[(df[0] == f"chr{chrom}") & (df[1] < end) & (df[2] > start)]

        # For each interval...
        for _, row in df.iterrows():

            # Get genomic feature
            feat = GenomicFeature(
                chrom=row[0],
                start=row[1],
                end=row[2],
                score=weight,
                feat_type=feat_type,
                feat_id=f"{row[0]}:{row[1]+1}-{row[2]}",
            )
            # evidence_feats[evidence][-1].append(feat)
            evidence_feats[-1].append(feat)

    return evidence_feats


def features_to_profile(start, end, features, default_score=1.0):
    """Return an numpy array of scores from features in the given
    range.

    Code adapted from darenillas "PhastCons.py"
    Original function: "compute_conservation_profile"
    URL: https://github.com/wassermanlab/GUD/blob/master/lib/ORCA/analysis/PhastCons.py
    """

    # Initialize
    profile = np.zeros(end - start)

    # For each feature...
    for feat in features:

        # Get indices
        i = max([0, feat.start - start])
        j = feat.end - start

        # Set 
        if hasattr(feat, "score"):
            profile[i:j] = feat.score
        else:
            profile[i:j] = default_score

    return profile


def profile_to_regions(chrom, start, end, genome, profile, min_score,
                       min_length, exons=[], stitch=False, feat_type="Region",
                       feat_id_prefix="Region", source=None,
                       reset_count=False):
    """Compute regions from profile of at least min. length
    that score above min.

    """

    # Initialize
    regions = []
    stitched_regions = []
    fine_regions = []

    # Get consecutive bp scoring above threshold
    consecutive = _get_consecutive_numbers(np.where(profile >= min_score)[0])

    # For each set of consecutive bp...
    for i, arr in enumerate(consecutive):

        # Get genomic feature
        feat = GenomicFeature(
            chrom=chrom,
            start=start+arr[0],
            end=start+arr[-1]+1,
            score=np.mean(profile[arr]),
            feat_type=feat_type,
            feat_id=f"{feat_id_prefix}{i+1}",
            qualifiers={"source": source, "genome": genome},
            profile=profile[arr],
        )
        regions.append(feat)

    # Stitch regions into larger regions
    if stitch:

        # For each region...
        for i in range(len(regions) - 1):

            # Get indices of stitched region
            stitched_idxs = _stitch_regions(start, end, regions, [i], profile,
                                            exons, min_score)

            # Get genomic feature
            if len(stitched_idxs) > 1:

                id = f"{regions[stitched_idxs[0]].id}-{stitched_idxs[-1]+1}"

                feat = GenomicFeature(
                    chrom=chrom,
                    start=regions[stitched_idxs[0]].start,
                    end=regions[stitched_idxs[-1]].end,
                    feat_type=feat_type,
                    feat_id=id,
                    qualifiers={"source": source},
                )
                feat.profile = profile[feat.start-start:feat.end-start]
                feat.score = np.mean(feat.profile)
                stitched_regions.append(feat)

    # Get long enough non-overlapping regions
    for region in sorted(regions + stitched_regions,
                         key=lambda x: x.end - x.start + 1, reverse=True):

        # Too short!
        if region.end - region.start < min_length:
            continue

        overlap = False

        for fine_region in fine_regions:

            if region.overlaps(fine_region):
                overlap = True
                break

        if not overlap:
            fine_regions.append(region)

    fine_regions.sort(key=lambda x: x.start)

    if reset_count:
        counter = 0
        for fr in fine_regions:
            fr.id=f"{feat_id_prefix}{counter+1}"
            counter += 1
  
    return fine_regions


def _stitch_regions(start, end, regions, region_idxs, profile, exons=[],
                     min_score=0.6):
    """Stitch regions into larger regions and return those that still
    would score above min. By default, regions separated by coding
    exons will be combined.

    NOTE: If there were three regions (A, B and C) and both A and B
    and B to C could be stitched (but not A, B and C), the function
    would return either A and B or B and C, whichever would results
    in a arger region.

    Code adapted from darenillas "PhastCons.py"
    Original function: "_combine_conserved_regions_exluding_exons"
    URL: https://github.com/wassermanlab/GUD/blob/master/lib/ORCA/analysis/PhastCons.py

    """

    # Initialize
    first_region = regions[region_idxs[0]]
    last_region = regions[region_idxs[-1]]
    next_region = regions[region_idxs[-1] + 1]

    # Do not combine regions separated by coding exons
    if exons:

        # Region : ----last-----next----
        # exon 1 : -ex------------------ : exon in between, False
        # exon 2 : ------------------ex- : exon in between, False
        # exon 3 : ----------ex--------- : exon in between, True
        # exon 4 : ----ex--------------- : exon in between, False
        # exon 5 : ---------------ex---- : exon in between, False
        # exon 6 : -------ex------------ : exon in between, True
        # exon 7 : ------------ex------- : exon in between, True

        for exon in exons:

            # exons 2, 5 (i.e., no exon in between)
            if exon.start > next_region.start:
                continue

            # exons 1, 4 (i.e., no exon in between)
            if exon.end < last_region.end:
                continue

            # exons 3, 6, 7 (i.e., exon in between)
            if exon.end >= last_region.end and exon.start <= next_region.start:
                return region_idxs

    # Stitched regions must score above threshold
    s = np.mean(profile[first_region.start-start:next_region.end-start])
    if s < min_score:
        return region_idxs

    # Continue stitching regions
    region_idxs.append(region_idxs[-1] + 1)
    if region_idxs[-1] + 1 < len(regions):
        return _stitch_regions(start, end, regions, region_idxs, profile,
                               exons, min_score)

    return region_idxs


def _get_consecutive_numbers(arr):
    """
    Function to get sets of consecutive numbers in a NumPy array
    :param arr: 1D NumPy array
    :return: list of NumPy arrays
    """

    arr = np.split(arr, np.where(np.diff(arr) != 1)[0] + 1)

    return arr


def _get_interval_tfs(chrom, start, end, genome):
    """
    Function to get the JASPAR TFBSs of an interval from UCSC
    :param chrom: str, chromosome name
    :param start: int, start coordinate (0-based)
    :param end: int, end coordinate
    :param genome: str, genome assembly in UCSC Genome Browser
    :return: Bio Seq, sequence
    """

    # Get URL
    url = "https://api.genome.ucsc.edu/getData/" + \
          "track?track=jaspar2022;maxItemsOutput=-1;" + \
         f"genome={genome};chrom=chr{chrom};start={start};end={end}"

    # Get response
    response = requests.get(url)

    return set([tf.upper() for r in response.json()["jaspar2022"] \
                for tf in r["TFName"].split("::")])


def _liftover_regions(regions, from_genome, to_genome, dummy_dir="/tmp/"):
    """
    Function to get the JASPAR TFBSs of an interval from UCSC
    :param regions: list, regulatory regions as GenomicFeature
    :param from_genome: str, genome assembly in UCSC Genome Browser
    :param to_genome: str, genome assembly in UCSC Genome Browser
    :return: Bio Seq, sequence
    """

    # Initialize
    regexp = r"chr(\S{1,2})$"
    regions_lo = []
    base_name = os.path.basename(__file__)
    pid = os.getpid()

    # Get liftOver chain file
    chain_file = os.path.join(OnTargetUtils.liftover_dir,
        f"{from_genome}To{to_genome.capitalize()}.over.chain.gz")

    # BED file
    bed_file = os.path.join(dummy_dir, "%s.%s.bed" % (base_name, pid))
    handle = ParseUtils._get_file_handle(bed_file, "w")
    for r in regions:
        handle.write(f"chr{r.chrom}\t{r.start}\t{r.end}\t{r.id}\n")
    handle.close()

    # LiftOver
    minMatch = 0.1 if from_genome[:2] != to_genome[:2] else 0.95
    liftover_file = os.path.join(dummy_dir, "%s.%s.lo" % (base_name, pid))
    unmapped_file = os.path.join(dummy_dir, "%s.%s.unmap" % (base_name, pid))
    cmd = f"liftOver -minMatch={minMatch} " + \
          f"{bed_file} {chain_file} {liftover_file} {unmapped_file}"
    p = sp.Popen([cmd], stdout=sp.DEVNULL, stderr=sp.DEVNULL, shell=True)
    p.wait() # wait for child process to terminate

    # Get liftOver regions
    df = pd.read_table(liftover_file, header=None, index_col=3)
    for r in regions:
        if r.id in df.index:
            row = df.loc[[r.id]]
            m = re.search(regexp, row[0].to_string())
            feat = GenomicFeature(
                chrom=m.group(1),
                start=int(row[1]),
                end=int(row[2]),
                score=r.score,
                feat_type=r.type,
                feat_id=r.id,
                qualifiers=r.qualifiers,
                profile=r.profile,
            )
            coord = f"{from_genome}:chr{r.chrom}:{r.start+1}-{r.end}"
            feat.qualifiers.setdefault("original coordinate", coord)
            feat.qualifiers.setdefault("liftOver genome", to_genome)
            regions_lo.append(feat)

    # Remove
    os.remove(bed_file)
    os.remove(liftover_file)
    os.remove(unmapped_file)

    return regions_lo


if __name__ == "__main__":
    cli()