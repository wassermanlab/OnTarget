#!/usr/bin/env python2.7

# TODOs:
# 1) Implement FANTOM5 TSSs for promoter definition (line 521)

import os, sys, re
from array import array
import argparse
from binning import containing_bins, contained_bins
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
try: import ConfigParser
except: import configparser as ConfigParser
import json
from lxml import etree
import numpy
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
import shutil
import warnings

# Append GUD/OnTarget modules to path
module_path = os.path.join(os.path.dirname(__file__), os.pardir)
sys.path.append(os.path.join(module_path, os.pardir))
sys.path.append(os.path.join(module_path, "GUD"))

# Import from GUD/OnTarget
from GUD.ORM.gene import Gene
from GUD.ORM.chrom_size import ChromSize
from GUD.ORM.conservation import Conservation
from GUD.ORM.dna_accessibility import DnaAccessibility
from GUD.ORM.histone_modification import HistoneModification
from GUD.ORM.repeat_mask import RepeatMask
from GUD.ORM.tad import Tad
from GUD.ORM.tf_binding import TfBinding

from OnTarget import OTglobals
from OnTarget.lib.region import Region, compute_profile, compute_regions
from OnTarget.lib.transcript import Transcript

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(module_path, "config.ini")
config.read(config_file)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via command
    line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(description="describe what the script does...")

    # Optional args
    region_group = parser.add_mutually_exclusive_group(required=True)
    region_group.add_argument("--region", nargs=3, help="Genomic region (e.g. \"chr7 142720001 143120000\")", metavar=("CHROM", "START", "END"))
    region_group.add_argument("--region-file", help="BED file containing a list of genomic regions")

    parser.add_argument("-l", "--length", type=int, default=int(config.get("Parameters", "min_enhancer_length")),
        help="Min. enhancer length (default from \"config.ini\" = %s)" % config.get("Parameters", "min_enhancer_length"))

    feats = ["accessibility", "conservation", "enhancer", "histone", "tf"]
    feat_group = parser.add_mutually_exclusive_group()
    feat_group.add_argument("--include", default=[], nargs="*", choices=feats,
        help="ONLY include genomic features of that type (i.e. %s)" % ", ".join(feats), metavar="FEATURE")
    feat_group.add_argument("--exclude", default=[], nargs="*", choices=feats,
        help="Exclude genomic features of that type (i.e. %s)" % ", ".join(feats), metavar="FEATURE")

    parser.add_argument("--mask-exons", type=int, default=1, choices=[0, 1], help="Mask protein-coding exons (i.e. nucleotides overlapping exons are discarded; default = 1)")
    parser.add_argument("--mask-repeats", type=int, default=1, choices=[0, 1], help="Mask repeat regions (i.e. nucleotides overlapping repeats are discarded; default = 1)")

    sample_group = parser.add_mutually_exclusive_group()
    sample_group.add_argument("--sample", default=[], nargs="*", help="Sample (limits RR search to using features from cells or tissues of the given sample; e.g. \"brain\")")
    sample_group.add_argument("--sample-file", help="File containing a list of samples")

    parser.add_argument("--ubiquitous", type=int, default=1, choices=[0, 1], help="BE ubiquitous (i.e. use features from other samples; default = 1)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode (default = False)")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default=config.get("MySQL", "db"),
        help="Database name (e.g. \"mm10\"; default from \"config.ini\" = %s)" % config.get("MySQL", "db"))
    mysql_group.add_argument("-H", "--host", default=config.get("MySQL", "host"),
        help="Host name (e.g. \"ontarget.cmmt.ubc.ca\"; default from \"config.ini\" = %s)" % config.get("MySQL", "host"))
    mysql_group.add_argument("-P", "--port", default=config.get("MySQL", "port"),
        help="User name (e.g. \"5506\"; default from \"config.ini\" = %s)" % config.get("MySQL", "port"))
    mysql_group.add_argument("-u", "--user", default=config.get("MySQL", "user"),
        help="User name (e.g. \"ontarget_r\"; default from \"config.ini\" = %s)" % config.get("MySQL", "user"))

    # Output args
    out_group = parser.add_argument_group("output arguments")
    out_group.add_argument("-o", "--out-dir", default="./", help="Output directory (default = ./)")
    out_group.add_argument("--enhancers", type=int, default=1, choices=[0, 1], help="Output enhancers (in BED format; default = 1)")
    out_group.add_argument("--features", type=int, default=0, choices=[0, 1], help="Output genomic features (in BED format; default = 0)")
    out_group.add_argument("--matrix", type=int, default=0, choices=[0, 1], help="Output feature matrix (in CSV format; default = 0)")
    out_group.add_argument("--signal", type=int, default=0, choices=[0, 1], help="Output feature signal (in PNG format; default = 0)")

    # Handle errors
    args = parser.parse_args()

    if not args.enhancers and not args.features and not args.matrix:
        parser.error("missing arguments: specify at least one output (i.e. \"-e\", \"-f\" or \"-m\") or type \"-h\" for help")

    return args

def get_genomic_features(session, chrom, start, end, sample=[],
    bins=[], ubiquitous=True, verbose=False):
    """
    """

    # Initialize
    feats = {}
    feats.setdefault("exons", [])
    histone_types = config.get("Parameters", "histone_types").split(",")
#    repeat_classes = config.get("Parameters", "repeat_classes").split(",")

    # Get genes
    genes = Gene.select_by_bin_range(session, chrom,
        start, end, bins)
    # Add genes to feats
    feats.setdefault("genes", genes)
    # Verbose
    if verbose:
        print("\nGenes:")

    # Add exons to feats
    for g in genes:
        feats["exons"].extend(g.coding_exons)
        # Verbose
        if verbose:
            print(g)

    # Add sequence to feats
    feats.setdefault("sequence", get_sequence_from_ucsc(
        chrom, start, end, config.get("MySQL", "db")))

    # Get conservation
    conservation = Conservation.select_by_bin_range(session,
        chrom, start, end, bins)
    # Compute conservation profile
    conservation_profile = compute_profile(start, end,
        conservation)
    # Compute conserved regions
#    conserved_regions = compute_regions(chrom, start,
#        end, conservation_profile, exons=feats["exons"],
#        min_score=float(config.get("Parameters", "min_conservation")),
#        min_length=int(config.get("Parameters", "min_conservation_length")),
#        stitch=True, label="ConservedRegion", type="conserved region",
#        source=conservation[0].source_name)
    conserved_regions = compute_regions(chrom, start,
        end, conservation_profile, exons=[],
        min_score=float(config.get("Parameters", "min_conservation")),
        min_length=int(config.get("Parameters", "min_conservation_length")),
        stitch=True, label="ConservedRegion", type="conserved region",
        source=conservation[0].source_name)
#    # Add conservation profile to feats
#    feats.setdefault("conservation_profile",
#        conservation_profile)
    # Add conserved regions to feats
    feats.setdefault("conserved_regions",
        conserved_regions)
    # Verbose
    if verbose:
        print("\nConserved Regions:")
        for cr in conserved_regions: print(cr)

    # Get DNA accessibility
    dna_accessibility = DnaAccessibility.select_by_bin_range(
        session, chrom, start, end, sample=sample, bins=bins)
#    # If DNA accessibility feats could not be found...
#    if not dna_accessibility and sample and ubiquitous:
#        warnings.warn("\nCell-/tissue-specific DNA Accessibility features could not be found!\n\tUsing features from other cells/tissues instead...\n")
#        # ... Use any DNA accessibility feats
#        dna_accessibility = DnaAccessibility.select_by_bin_range(
#            session, chrom, start, end, bins=bins)
    # Add DNA accessibility to feats
    feats.setdefault("dna_accessibility",
        dna_accessibility)
    # Verbose
    if verbose:
        print("\nDNA Accessibility:")
        for da in dna_accessibility: print(da)

    # Get histone modifications
#    histone_modification = HistoneModification.select_by_bin_range(
#        session, chrom, start, end, histone_types=histone_types,
#        sample=sample, bins=bins)
    histone_modification = HistoneModification.select_by_bin_range(
        session, chrom, start, end, histone_types=[], sample=sample, bins=bins)
#    # If histone modification feats could not be found...
#    if not histone_modification and sample and ubiquitous:
#        warnings.warn("\nCell-/tissue-specific Histone Modification features could not be found!\n\tUsing features from other cells/tissues instead...\n")
#        # ... Use any histone modification feats
#        histone_modification = HistoneModification.select_by_bin_range(
#            session, chrom, start, end, histone_types=histone_types,
#            bins=bins)
    # Add histone modifications to feats
    feats.setdefault("histone_modification",
        histone_modification)
    # Verbose
    if verbose:
        print("\nHistone Modifications:")
        for hm in histone_modification: print(hm)

#    # Get repeat regions
#    repeats = RepeatMask.select_by_bin_range(session, chrom, start,
#        end, repeat_classes=repeat_classes, bins=bins)
#    # Add repeats to feats
#    feats.setdefault("repeats", repeats)
#    # Verbose
#    if verbose:
#        print("\nRepeat Regions:")
#        for r in repeats: print(r)

#    # Get TADs encompassing the gene
#    tads = Tad.select_encompasing_range(session, chrom, start, end,
#        sample=sample)
#    # If TADs could not be found...
#    if not tads:
#        # ... Use TADs overlapping the gene
#        tads = Tad.select_overlapping_range(session, chrom, start,
#            end, sample=sample)
#    # If TADs could not be found...
#    if not tads and sample and ubiquitous:
#        warnings.warn("\nCell-/tissue-specific TADs could not be found!\n\tUsing TADs from other cells/tissues instead...\n")
#        # ... Use any TADs encompassing the gene
#        tads = Tad.select_encompasing_range(session, chrom,
#            start, end)
#        # If TADs could not be found...
#        if not tads:
#            # ... Use any TADs overlapping the gene
#            tads = Tad.select_overlapping_range(session, chrom,
#                start, end)
#    # Add TADs to feats
#    feats.setdefault("tads", tads)
#    # Verbose
#    if verbose:
#        print("\nTADs:")
#        for tad in tads: print(tad)

    # Get TF-binding (use ALL ReMap)
    tf_binding = TfBinding.select_by_bin_range(
        session, chrom, start, end, bins=bins)
#    # If TF-binding feats could not be found...
#    if not tf_binding:
#        warnings.warn("\nCell-/tissue-specific TF-Binding features could not be found!\n\tUsing features from other cells/tissues instead...\n")
#        # ... Use any TF-binding feats
#        tf_binding = TfBinding.select_by_bin_range(session,
#            chrom, start, end, bins=bins)
    # Add TF-binding to feats
    feats.setdefault("tf_binding", tf_binding)
    # Verbose
    if verbose:
        print("\nTF-Binding:")
        for tf in tf_binding: print(tf)

    return feats

def get_sequence_from_ucsc(chrom=None, start=None, end=None, genome=None):
    """
    Retrieve a DNA sequence from UCSC.

    Note: the function assumes that the start position
    is 1-based.
    """

    # Initialize
    sequence = ""
    url = "http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=%s:%s,%s" % (
        genome, chrom, start, end)

    # Get XML
    xml = etree.parse(url, parser=etree.XMLParser())
    # Get sequence
    sequence = xml.xpath("SEQUENCE/DNA/text()")[0].replace("\n", "")

    return Seq(sequence, unambiguous_dna)

def filter_genomic_features(feats, include=[], exclude=[]):
    """
    """

    # For each feat...
    for feat in frozenset(feats.keys()):
        # Exclude conserved regions
        if feat == "conserved_regions":
            if include:
                if "conservation" not in include: del feats[feat]
            elif "conservation" in exclude: del feats[feat]
        # Exclude DNA accessibility
        elif feat == "dna_accessibility":
            if include:
                if "accessibility" not in include: del feats[feat]
            elif "accessibility" in exclude: del feats[feat]
        # Exclude histone modifications
        elif feat == "histone_modification":
            if include:
                if "histone" not in include: del feats[feat]
            elif "histone" in exclude: del feats[feat]
        # Exclude histone modifications
        elif feat == "tf_binding":
            if include:
                if "tf" not in include: del feats[feat]
            elif "tf" in exclude: del feats[feat]

    return feats

def get_evidence_and_weights(feats):
    """
    """

    # Initialize
    evidence = {}
    weights = {}

#    # Get conservation profile evidence
#    if "conservation_profile" in feats:
#        evidence.setdefault("conservation_profile",
#            feats["conservation_profile"])
#        # If weight is defined for conservation profile...
#        if config.has_option("Weights", "conservation_profile"):
#            # Assign defined weight
#            weights.setdefault("conservation_profile",
#                float(config.get("Weights", "conservation_profile")))
#        # ... Else...
#        else:
#            warnings.warn("\nWeight for \"conservation_profile\" not defined (refer to \"config.ini\" >>> \"Weights\")!\n\tIgnoring evidence...\n")
#            # Assign a weight of 0
#            weights.setdefault("conservation_profile", 0.0)

    # Get conserved regions evidence
    if "conserved_regions" in feats:
        evidence.setdefault("conserved_regions",
            feats["conserved_regions"])
#        # If weight is defined for conserved regions...
#        if config.has_option("Weights", "conserved_regions"):
#            # Assign defined weight
#            weights.setdefault("conserved_regions",
#                float(config.get("Weights", "conserved_regions")))
#        else:
#            warnings.warn("\nWeight for \"conserved_regions\" not defined (refer to \"config.ini\" >>> \"Weights\")!\n\tIgnoring evidence...\n")
#            # Assign a weight of 0
#            weights.setdefault("conserved_regions", 0.0)

    # For each DNA accessibility feat...
    if "dna_accessibility" in feats:
        for feat in feats["dna_accessibility"]:
            # Unwind feats by experiment
            evidence.setdefault(feat.experiment_type, [])
            evidence[feat.experiment_type].append(feat)
#            # If weight is defined for experiment...
#            if config.has_option("Weights.DnaAccessibility", feat.experiment_type):
#                # Assign defined weight
#                weights.setdefault(feat.experiment_type,
#                    float(config.get("Weights.DnaAccessibility", feat.experiment_type)))
#            # ... Instead, if general weight is defined for DNA accessibility...
#            elif config.has_option("Weights", "dna_accessibility"):
#                warnings.warn("\nWeight for \"%s\" not defined (refer to \"config.ini\" >>> \"Weights.DnaAccessibility\")!\n\tUsing general \"dna_accessibility\" weight instead...\n" % 
#                    feat.experiment_type)
#                # Assign defined weight
#                weights.setdefault(feat.experiment_type,
#                    float(config.get("Weights", "dna_accessibility")))
#            # ... Else...
#            else:
#                warnings.warn("\nWeight for \"dna_accessibility\" not defined (refer to \"config.ini\" >>> \"Weights\")!\n\tIgnoring evidence...\n")
#                # Assign a weight of 0
#                weights.setdefault(feat.experiment_type, 0.0)

    # For each histone modification feat...
    if "histone_modification" in feats:
        for feat in feats["histone_modification"]:
            # Unwind feats by histone
            evidence.setdefault(feat.histone_type, [])
            evidence[feat.histone_type].append(feat)
#            # If weight is defined for histone...
#            if config.has_option("Weights.HistoneModification", feat.histone_type):
#                # Assign defined weight
#                weights.setdefault(feat.histone_type,
#                    float(config.get("Weights.HistoneModification", feat.histone_type)))
#            # ... Instead, if general weight is defined for histone modifications...
#            elif config.has_option("Weights", "histone_modification"):
#                warnings.warn("\nWeight for \"%s\" not defined (refer to \"config.ini\" >>> \"Weights.HistoneModification\")!\n\tUsing general \"histone_modification\" weight instead...\n" % 
#                    feat.histone_type)
#                # Assign defined weight
#                weights.setdefault(feat.histone_type,
#                    float(config.get("Weights", "histone_modification")))
#            # ... Else...
#            else:
#                warnings.warn("\nWeight for \"histone_modification\" not defined (refer to \"config.ini\" >>> \"Weights\")!\n\tIgnoring evidence...\n")
#                # Assign a weight of 0
#                weights.setdefault(feat.histone_type, 0.0)

    # For each TF binding feat...
    if "tf_binding" in feats:
        for feat in feats["tf_binding"]:
            # Unwind feats by TF
            evidence.setdefault(feat.tf_name, [])
            evidence[feat.tf_name].append(feat)
#            # If weight is defined for TF...
#            if config.has_option("Weights.TfBinding", feat.tf_name):
#                # Assign defined weight
#                weights.setdefault(feat.tf_name,
#                    float(config.get("Weights.TfBinding", feat.tf_name)))
#            # ... Instead, if general weight is defined for TF-binding...
#            elif config.has_option("Weights", "tf_binding"):
#                warnings.warn("\nWeight for \"%s\" not defined (refer to \"config.ini\" >>> \"Weights.TfBinding\")!\n\tUsing general \"tf_binding\" weight instead...\n" %
#                    feat.tf_name)
#                # Assign defined weight
#                weights.setdefault(feat.tf_name,
#                    float(config.get("Weights", "tf_binding")))
#            # ... Else...
#            else:
#                warnings.warn("\nWeight for \"tf_binding\" not defined (refer to \"config.ini\" >>> \"Weights\")!\n\tIgnoring evidence...\n")
#                # Assign a weight of 0
#                weights.setdefault(feat.histone_type, 0.0)

    return evidence, weights

def compute_feature_matrix(start, end, evidence, weights):
    """
    Build a matrix of scores where each row represents a type
    of evidence and each column a genomic position. For each
    position, for each type of evidence, if \"discrete\", set
    the score at that position to 1 and, if \"continuous\", set
    the score at that position as defined in the evidence array.
    Scores are multiplied by the corresponding evidence weight.
    Background is set to 0.

    Weigths and evidence are dictionaries where the keys are
    type of evidence.
    
    NOTE: keys must agree between the weights and evidence!
    """

    # Initialize
    feat_matrix = None
    matrix_rows = []

    # For each feature type (e.g. DNase-seq, H3K4me3, etc.)...
    for feat_type in evidence:
        # If evidence is continuous...
        if isinstance(evidence[feat_type], array):
            # Set the feature vector to the evidence array
            feat_vector = numpy.array(evidence[feat_type])
#            # Multiply the feature vector by the weight
#            feat_vector *= weights[feat_type]
        # ... Instead, if evidence is discrete... #
        else:
            # Initialize the feature vector with 0s
            feat_vector = numpy.zeros(end - start, numpy.float)
            # For each feature...
            for feat in evidence[feat_type]:
                # Set each vector position to the feature
                # type weight
#                feat_vector[feat.start - start:
#                    feat.end - start] = weights[feat_type]
                feat_vector[feat.start - start:
                    feat.end - start] = 1.0

        # If first feature vector...
        if feat_matrix is None:
            # Initialize the matrix with the vector
            feat_matrix = feat_vector
        # ... Else...
        else:
            # Add vector to the matrix
            feat_matrix = numpy.vstack((feat_matrix, feat_vector))
        # Add feature type to rows
        matrix_rows.append(feat_type)

    return feat_matrix, matrix_rows

def compute_enhancers(chrom, start, end, feats, weights,
    feat_matrix, mask_exons=True, mask_repeats=True, verbose=False):
    """
    """

    # Initialize
    min_enhancer_len = int(
        config.get("Parameters", "min_enhancer_length"))
    perc_score_thresh = float(
        config.get("Parameters", "perc_score_thresh"))

    # Get the max. possible score of any enhancer
    max_score = sum(filter(lambda x: x > 0, weights.values()))

    # Summarize the matrix into a vector
    score_vector = numpy.sum(feat_matrix, axis=0)

    # Mask out protein coding exons and repeat regions
    if mask_exons:
        score_vector *= create_feature_mask(start, end, feats["exons"])
    if mask_repeats:
        score_vector *= create_feature_mask(start, end, feats["repeats"])

    # Normalize per-nucleotide score
    score_vector /= max_score

    # Use the top score percentile to determine a score
    # threshold used to define enhancer boundaries
    sorted_vector = numpy.sort(score_vector)
    top_perc_pos = sorted_vector.size * perc_score_thresh
    score_thresh = sorted_vector[int(sorted_vector.size - top_perc_pos)]

    # Compute conserved regions
    enhancers = compute_regions(chrom, start,
        end, array("f", score_vector), exons=feats["exons"],
        min_score=score_thresh, min_length=min_enhancer_len,
        label="RegulatoryRegion", type="Enhancer",
        source="OnTarget")

    # For each enhancer...
    for enhancer in enhancers:
        # Re-scale score between 0 and 1000
        # for display in the UCSC Genome Browser
        enhancer.score *= 1000
        # For each gene...
        for gene in feats["genes"]:
            # Skip if enhancer does not overlap TSS
            if enhancer.end < gene.txStart: continue
            if enhancer.start > gene.txStart: continue

            ##################
            # FANTOM5 TSS!!! #
            ##################

            # Identify promoters
            enhancer.type = "%s.Promoter" % gene.name2
            break

    return enhancers

def create_feature_mask(start, end, feats):

    # Initialize
    mask = numpy.ones(end - start, numpy.int)

    # For each feature...
    for feat in feats:
        # If feature overlaps the region...
        if feat.end > start and feat.start < end:
            mask_start = feat.start - start
            mask_end   = feat.end - start
            # If feat overlaps the region partially...
            if mask_start < 0:
                mask_start = 0
            if mask_end > end - start:
                mask_end = end - start
            # Set feat positions within vector to 0
            mask[mask_start : mask_end] = 0

    return mask

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Fetch genomic regions
    if args.region:
        regions = [[args.region[0], int(args.region[1]), int(args.region[2])]]
    elif args.region_file:
        region_file = os.path.abspath(args.region_file)
        regions = [[r[0], int(r[1]), int(r[2])] for r in OTglobals.parse_tsv_file(region_file)]

    # Fetch samples
    if args.sample:
        samples = args.sample
    elif args.sample_file:
        sample_file = os.path.abspath(args.sample_file)
        samples = [s for s in OTglobals.parse_file(sample_file)]
    else: samples = []

    # Establish a MySQL session
    db_name = "mysql://{}:@{}:{}/{}".format(args.user, args.host, args.port, args.db)
    try:
        engine = create_engine(db_name, echo=False)
        session = Session(engine)
    except:
        raise ValueError("Cannot connect to MySQL: %s" % db_name)

    # Get chromosome sizes
    chrom_sizes = ChromSize.chrom_sizes(session)

    # For each region... #
    for chrom, start, end in regions:
        try:
            # Validate region
            if chrom not in chrom_sizes:
                raise ValueError("\"{}\" is not a valid chrom!".format(chrom))
            if start <= 0:
                raise ValueError("\"{}\" is not a valid start position! Use \"1\" or higher.".format(start))
            if end > chrom_sizes[chrom]:
                raise ValueError("\"{}\" is not a valid end position! Use \"{}\" or lower.".format(end, chrom_sizes[chrom]))
            if end - start + 1 < args.length:
                raise ValueError("Insuficient region length! Min. region length is \"{}\".".format(args.length))

            # Zero the region start
            zeroed_start = start - 1

            # Pre-compute bins to speed up database selections
            # instead of computing them every time
            bins = list(set(containing_bins(zeroed_start, end) +\
                       contained_bins(zeroed_start, end)))

            # Get genomic features
            feats = get_genomic_features(session, chrom, zeroed_start,
                end, samples, bins, args.ubiquitous, args.verbose)

    #        # Filter genomic features
    #        feats = filter_genomic_features(feats, args.include,
    #            args.exclude)

            # Get evidence
            evidence, weights = get_evidence_and_weights(feats)

            # Compute matrix
            feat_matrix, matrix_rows = compute_feature_matrix(
                zeroed_start, end, evidence, weights)

            # Compute enhancers
            enhancers = compute_enhancers(chrom, zeroed_start, end,
                feats, weights, feat_matrix, args.mask_exons,
                args.mask_repeats, args.verbose)

            # Create output dir
            out_dir = os.path.abspath(args.out_dir)
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
    
            # Print enhancers
            if args.enhancers:
                sequence = str(feats["sequence"])
                enhancers_bed = os.path.join(out_dir, "%s:%s-%s.enhancers.bed" %
                    (chrom, start, end))
                enhancers_fasta = os.path.join(out_dir, "%s:%s-%s.enhancers.fa" %
                    (chrom, start, end))
                for e in enhancers:
                    OTglobals.write(enhancers_bed, e)
                    sub_sequence = sequence[e.start_1_based - zeroed_start:e.end - zeroed_start + 1]
                    OTglobals.write(enhancers_fasta, ">%s:%s-%s\n%s" % (chrom,
                        e.start_1_based, e.end, sub_sequence))

    #        # Print features
    #        if args.features:
    #            # For each feature...
    #            for feat_type in feats:
    #                # BED files
    #                if feat_type in ["conserved_regions", "dna_accessibility", "exons", "genes", "histone_modification", "repeats", "tads", "tf_binding"]:
    #                    feat_file = os.path.join(out_dir, "%s:%s-%s.%s.bed" %
    #                        (chrom, start, end, feat_type))
    #                    for feat in feats[feat_type]:
    #                        OTglobals.write(feat_file, feat)
    #                # FASTA file
    #                if feat_type == "sequence":
    #                    feat_file = os.path.join(out_dir, "%s:%s-%s.fa" %
    #                        (chrom, start, end))
    #                    OTglobals.write(feat_file, ">%s:%s-%s\n%s" % (chrom, start, end, feats["sequence"]))

            # Print matrix
            if args.matrix:
                if len(matrix_rows) == 1: continue
                # Initalize
                matrix = {}
                matrix_file = os.path.join(out_dir, "%s:%s-%s.matrix.csv" %
                    (chrom, start, end))
                # For each row...
                for i in range(len(matrix_rows)):
                    OTglobals.write(matrix_file, "%s,%s" % (matrix_rows[i], ",".join(
                        map(str, feat_matrix[i].tolist()))))

    #        # Print matrix
    #        if args.signal:
    #            # Initalize
    #            matrix = {}
    #            matrix_file = os.path.join(out_dir, "%s:%s-%s.matrix.csv" %
    #                (chrom, start, end))
    #            # For each row...
    #            for i in range(len(matrix_rows)):
    #                OTglobals.write(matrix_file, "%s,%s" % (matrix_rows[i], ",".join(
    #                    map(str, feat_matrix[i].tolist()))))
        except:
            print("Skipping region: %s:%s-%s" % (chrom, start, end))