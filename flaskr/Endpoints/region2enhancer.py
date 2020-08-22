import requests
from lxml import etree
import os, sys, re
from array import array
from binning import containing_bins, contained_bins
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import json
import numpy
import shutil
import warnings


# region should be an array of chromesome, start, stop
def region2enhancer(remote_host, database, region, length, feats, mask_exons, mask_repeats, ubiquitous, samples=[]):
    URL = remote_host + '/api/v1/' + database
    feats = ["accessibility", "conservation", "enhancer", "histone", "tf"]
    # Fetch genomic regions
    regions = [region[0], int(region[1]), int(region[2])]

    # Get chromosome sizes
    chromURL = URL+ '/chroms'
    try:
        chromData = requests.get(url=chromURL).json()
    except:
        print("error connecting to GUD")
    chrom_sizes = {}
    for item in chromData['results']:
        chrom_sizes[item['chrom']] = item['size']
    # For each region... #
    for chrom, start, end in regions:
        try:
            # Validate region
            if chrom not in chrom_sizes:
                raise ValueError("\"{}\" is not a valid chrom!".format(chrom))
            if start <= 0:
                raise ValueError("\"{}\" is not a valid start position! Use \"1\" or higher.".format(start))
            if end > chrom_sizes[chrom]:
                raise ValueError(
                    "\"{}\" is not a valid end position! Use \"{}\" or lower.".format(end, chrom_sizes[chrom]))
            if end - start + 1 < length:
                raise ValueError("Insuficient region length! Min. region length is \"{}\".".format(length))

            # Zero the region start
            zeroed_start = start - 1

            # Pre-compute bins to speed up database selections
            # instead of computing them every time
            bins = list(set(containing_bins(zeroed_start, end) + \
                            contained_bins(zeroed_start, end)))

            # Get genomic features
            feats = get_genomic_features(URL, chrom, zeroed_start,
                                         end, samples, bins, ubiquitous, verbose)

            # Get evidence
            evidence, weights = get_evidence_and_weights(feats)

            # Compute matrix
            feat_matrix, matrix_rows = compute_feature_matrix(
                zeroed_start, end, evidence, weights)

            # Compute enhancers
            enhancers = compute_enhancers(chrom, zeroed_start, end,
                                          feats, weights, feat_matrix, mask_exons,
                                          mask_repeats, verbose)

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

        except:
            print("Skipping region: %s:%s-%s" % (chrom, start, end))


def get_genomic_features(URL, chrom, start, end, sample=[],
                         bins=[], ubiquitous=True, verbose=False):
    """
    """

    # Initialize
    feats = {}
    feats.setdefault("exons", [])
    histone_types = config.get("Parameters", "histone_types").split(",")
    #    repeat_classes = config.get("Parameters", "repeat_classes").split(",")

    # Get genes
    genes=requests.get(url=URL+"/genes", params={'chrom':chrom,'start':start,'end':end,'location':'overlapping'}).json()
    if genes["results"]:
        genes=genes["results"]
    else:
        raise ValueError("cannot connect to DB")
    # Add genes to feats
    feats.setdefault("genes", genes)
    # Verbose
    if verbose:
        print("\nGenes:")

    # Add exons to feats
    for g in genes:
        coding_exons=zip(g['qualifiers']['exon_starts'],g['qualifiers']['exon_end'])
        feats["exons"].extend(coding_exons)
        # Verbose
        if verbose:
            print(g)

    # Add sequence to feats
    feats.setdefault("sequence", get_sequence_from_ucsc(
        chrom, start, end, config.get("MySQL", "db")))

    # Get conservation
    conservation=requests.get(url=URL+"/conservation", params={'chrom':chrom,'start':start,'end':end,'location':'overlapping'}).json()
    if conservation["results"]:
        conservation=conservation["results"]
    else:
        raise ValueError("cannot connect to DB")
    # Compute conservation profile
    conservation_profile = compute_profile(start, end,
                                           conservation)
    # Compute conserved regions
    conserved_regions = compute_regions(chrom, start,
                                        end, conservation_profile, exons=[],
                                        min_score=float(config.get("Parameters", "min_conservation")),
                                        min_length=int(config.get("Parameters", "min_conservation_length")),
                                        stitch=True, label="ConservedRegion", type="conserved region",
                                        source=conservation[0].source_name)
    # Add conserved regions to feats
    feats.setdefault("conserved_regions",
                     conserved_regions)
    # Verbose
    if verbose:
        print("\nConserved Regions:")
        for cr in conserved_regions: print(cr)

    # Get DNA accessibility
    dna_accessibility=requests.get(url=URL+"/dna_accessibility", params={'chrom':chrom,'start':start,'end':end,'location':'overlapping'}).json()
    if dna_accessibility["results"]:
        dna_accessibility=dna_accessibility["results"]
    else:
        raise ValueError("cannot connect to DB")
    # Add DNA accessibility to feats
    feats.setdefault("dna_accessibility",
                     dna_accessibility)
    # Verbose
    if verbose:
        print("\nDNA Accessibility:")
        for da in dna_accessibility: print(da)

    # Get histone modifications
    # histone_modification = HistoneModification.select_by_bin_range(
    #     session, chrom, start, end, histone_types=[], sample=sample, bins=bins)
    histone_modification=requests.get(url=URL+"/histone_modifications", params={'chrom':chrom,'start':start,'end':end,'histone_types':histone_types,'sample':sample,'location':'overlapping'}).json()
    if histone_modification["results"]:
        histone_modification=histone_modification["results"]
    else:
        raise ValueError("cannot connect to DB")
    # Add histone modifications to feats
    feats.setdefault("histone_modification",
                     histone_modification)
    # Verbose
    if verbose:
        print("\nHistone Modifications:")
        for hm in histone_modification: print(hm)

    # Get TF-binding (use ALL ReMap)
    tf_binding=requests.get(url=URL+"/tf_binding", params={'chrom':chrom,'start':start,'end':end,'location':'overlapping'}).json()
    if tf_binding["results"]:
        tf_binding=tf_binding["results"]
    else:
        raise ValueError("cannot connect to DB")
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


def get_evidence_and_weights(feats):
    """
    """

    # Initialize
    evidence = {}
    weights = {}

    # Get conserved regions evidence
    if "conserved_regions" in feats:
        evidence.setdefault("conserved_regions",
                            feats["conserved_regions"])

    # For each DNA accessibility feat...
    if "dna_accessibility" in feats:
        for feat in feats["dna_accessibility"]:
            # Unwind feats by experiment
            evidence.setdefault(feat.experiment_type, [])
            evidence[feat.experiment_type].append(feat)
    # For each histone modification feat...
    if "histone_modification" in feats:
        for feat in feats["histone_modification"]:
            # Unwind feats by histone
            evidence.setdefault(feat.histone_type, [])
            evidence[feat.histone_type].append(feat)
    # For each TF binding feat...
    if "tf_binding" in feats:
        for feat in feats["tf_binding"]:
            # Unwind feats by TF
            evidence.setdefault(feat.tf_name, [])
            evidence[feat.tf_name].append(feat)

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
            mask_end = feat.end - start
            # If feat overlaps the region partially...
            if mask_start < 0:
                mask_start = 0
            if mask_end > end - start:
                mask_end = end - start
            # Set feat positions within vector to 0
            mask[mask_start: mask_end] = 0

    return mask
