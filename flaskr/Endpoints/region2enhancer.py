import requests
from lxml import etree
import os, sys, re
from array import array
from binning import containing_bins, contained_bins
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
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
    chromURL = URL + '/chroms'
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
    genes = requests.get(url=URL + "/genes",
                         params={'chrom': chrom, 'start': start, 'end': end, 'location': 'overlapping'}).json()
    if genes["results"]:
        genes = genes["results"]
    else:
        raise ValueError("cannot connect to DB")
    # Add genes to feats
    feats.setdefault("genes", genes)
    # Verbose
    if verbose:
        print("\nGenes:")

    # Add exons to feats
    for g in genes:
        coding_exons = zip(g['qualifiers']['exon_starts'], g['qualifiers']['exon_end'])
        feats["exons"].extend(coding_exons)
        # Verbose
        if verbose:
            print(g)

    # Add sequence to feats
    feats.setdefault("sequence", get_sequence_from_ucsc(
        chrom, start, end, config.get("MySQL", "db")))

    # Get conservation
    conservation = requests.get(url=URL + "/conservation",
                                params={'chrom': chrom, 'start': start, 'end': end, 'location': 'overlapping'}).json()
    if conservation["results"]:
        conservation = conservation["results"]
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
    dna_accessibility = requests.get(url=URL + "/dna_accessibility", params={'chrom': chrom, 'start': start, 'end': end,
                                                                             'location': 'overlapping'}).json()
    if dna_accessibility["results"]:
        dna_accessibility = dna_accessibility["results"]
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
    histone_modification = requests.get(url=URL + "/histone_modifications",
                                        params={'chrom': chrom, 'start': start, 'end': end,
                                                'histone_types': histone_types, 'sample': sample,
                                                'location': 'overlapping'}).json()
    if histone_modification["results"]:
        histone_modification = histone_modification["results"]
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
    tf_binding = requests.get(url=URL + "/tf_binding",
                              params={'chrom': chrom, 'start': start, 'end': end, 'location': 'overlapping'}).json()
    if tf_binding["results"]:
        tf_binding = tf_binding["results"]
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


def compute_profile(start, end, features):
    """
    Return a profile (i.e. array of scores) from features
    in the given region (i.e. start, end).

    Code adapted from darenillas "PhastCons.py"
    Original function: "compute_conservation_profile"
    URL: https://github.com/wassermanlab/GUD/blob/master/lib/ORCA/analysis/PhastCons.py
    """

    # Initialize
    scores = {}
    profile = array("f")

    # Extract per-position scores
    for feat in features:
        for i in range(feat.start, feat.end):
            if hasattr(feat, "score"):
                scores.setdefault(i, feat.score)
            # Assign position a score of 1.0
            else:
                scores.setdefault(i, 1.0)

    # Fill gaps with 0s
    for i in range(start, end):
        if i in scores:
            profile.append(scores[i])
        else:
            profile.append(0.0)

    return profile


def compute_regions(chrom, start, end, profile, exons=[],
                    min_score=0.6, min_length=10, stitch=False,
                    label="Region", type="region", source=None):
    """
    Compute regions from profile of at least min. length
    that score above min.
    """

    counter = 0
    start_idx = 0
    end_idx = end - start
    region_start_idx = None
    region_end_idx = None
    regions = []
    stitched_regions = []
    fine_regions = []
    # This should call 6 conserved regions, and merge 2-3-4, and 5-6
    # profile = array('f', [6.80556863e-01, 9.43449691e-01, 8.37105628e-01, 6.53790967e-01,
    #                            7.94108937e-01, 8.92694225e-01, 8.70977562e-01, 6.24933764e-01,
    #                            9.30884657e-01, 7.32830332e-01, 7.34998617e-01, 6.10970434e-01,
    #                            9.99887508e-01, 9.06141433e-01, 8.07754123e-01, 7.44507877e-01,
    #                            6.20665811e-01, 7.91991350e-01, 6.02578636e-01, 8.60206886e-01,
    #                            7.38006268e-01, 8.99331196e-01, 6.37047916e-01, 6.39402431e-01,
    #                            7.65803396e-01, 1.55889807e-01, 2.21438493e-02, 3.04482285e-02,
    #                            3.19881397e-04, 1.42229888e-01, 2.93923127e-02, 8.61047507e-02,
    #                            1.04912087e-01, 1.35948892e-01, 9.54951895e-04, 1.75901521e-01,
    #                            2.85730217e-02, 5.73076556e-03, 1.10217622e-01, 2.35927857e-01,
    #                            2.40422533e-01, 1.82789568e-01, 1.88084750e-01, 1.30239308e-01,
    #                            1.22399532e-01, 2.27598554e-01, 8.87512698e-02, 1.45089697e-01,
    #                            2.31214382e-01, 2.19829281e-01, 2.20741445e-01, 1.01574773e-01,
    #                            3.00931665e-02, 4.15615321e-02, 4.61427200e-03, 2.21379064e-01,
    #                            2.88507097e-02, 8.52727889e-02, 2.03898513e-01, 2.28187987e-01,
    #                            6.27985068e-01, 8.04470639e-01, 7.22583865e-01, 9.11400757e-01,
    #                            6.79460397e-01, 6.38220219e-01, 8.12702963e-01, 7.62883328e-01,
    #                            6.44603354e-01, 8.43406180e-01, 6.17595060e-01, 7.05057351e-01,
    #                            8.56613324e-01, 6.20265629e-01, 6.64293609e-01, 4.98575204e-02,
    #                            2.86441430e-01, 1.06564167e-01, 2.78406493e-01, 5.60276294e-02,
    #                            6.24925026e-01, 6.78700171e-01, 7.43820261e-01, 7.15024611e-01,
    #                            8.05726999e-01, 8.88167621e-01, 6.88128705e-01, 6.59038054e-01,
    #                            8.89519960e-01, 9.88556404e-01, 7.91321195e-01, 8.90755524e-01,
    #                            9.14818075e-01, 8.34835039e-01, 6.40025535e-01,
    #                            1.22399532e-01, 2.27598554e-01, 8.87512698e-02, 1.45089697e-01,
    #                            7.94108937e-01, 8.92694225e-01, 8.70977562e-01, 6.24933764e-01,
    #                            9.30884657e-01, 7.32830332e-01, 7.34998617e-01, 6.10970434e-01,
    #                            9.99887508e-01, 9.06141433e-01, 8.07754123e-01, 7.44507877e-01,
    #                            6.20665811e-01, 7.91991350e-01, 6.02578636e-01, 8.60206886e-01,
    #                            7.38006268e-01, 8.99331196e-01, 6.37047916e-01, 6.39402431e-01,
    #                            7.65803396e-01, 1.55889807e-01, 2.21438493e-02, 3.04482285e-02,
    #                            3.19881397e-04, 1.42229888e-01, 2.93923127e-02, 8.61047507e-02,
    #                            1.04912087e-01, 1.35948892e-01, 9.54951895e-04, 1.75901521e-01,
    #                            2.85730217e-02, 5.73076556e-03, 1.10217622e-01, 2.35927857e-01,
    #                            2.40422533e-01, 1.82789568e-01, 1.88084750e-01, 1.30239308e-01,
    #                            1.22399532e-01, 2.27598554e-01, 8.87512698e-02, 1.45089697e-01,
    #                            2.31214382e-01, 2.19829281e-01, 2.20741445e-01, 1.01574773e-01,
    #                            3.00931665e-02, 4.15615321e-02, 4.61427200e-03, 2.21379064e-01,
    #                            2.88507097e-02, 8.52727889e-02, 2.03898513e-01, 2.28187987e-01,
    #                            6.27985068e-01, 8.04470639e-01, 7.22583865e-01, 9.11400757e-01,
    #                            6.79460397e-01, 6.38220219e-01, 8.12702963e-01, 7.62883328e-01,
    #                            6.44603354e-01, 8.43406180e-01, 6.17595060e-01, 7.05057351e-01,
    #                            8.56613324e-01, 6.20265629e-01, 6.64293609e-01, 4.98575204e-02,
    #                            2.86441430e-01, 1.06564167e-01, 2.78406493e-01, 5.60276294e-02,
    #                            6.24925026e-01, 6.78700171e-01, 7.43820261e-01, 7.15024611e-01,
    #                            8.05726999e-01, 8.88167621e-01, 6.88128705e-01, 6.59038054e-01,
    #                            8.89519960e-01, 9.88556404e-01, 7.91321195e-01, 8.90755524e-01,
    # 9.14818075e-01, 8.34835039e-01, 6.40025535e-01])
    #    end_idx = len(profile) - 1
    #    start = start_idx + 1
    #    end = end_idx + 1
    #    exons = []

    # For each exon...
    for exon in exons:
        # Initialize mask
        mask_start_idx = exon.start - start
        mask_end_idx = exon.end - start
        # If exon does not overlap region...
        if mask_end_idx <= start_idx or mask_start_idx >= end_idx:
            continue
        # ... Instead, if exon overlaps region...
        if mask_start_idx < start_idx:
            mask_start_idx = start_idx
        if mask_end_idx > end_idx:
            mask_end_idx = end_idx
        # Mask exons
        profile[mask_start_idx:mask_end_idx] = array("f",
                                                     [0.0] * (mask_end_idx - mask_start_idx))

    # For each position...
    for i in range(len(profile)):
        # If position scores above threshold...
        if profile[i] >= min_score:
            # Initialize region end index
            region_end_idx = i
            # Initialize region start index
            if region_start_idx is None:
                region_start_idx = region_end_idx
        # ... Instead, if position scores below threshold...
        else:
            # If region start/end indices are defined...
            if region_start_idx is not None:
                # Initialize region coordinates
                region_start = start + region_start_idx
                region_end = start + region_end_idx
                # Get region profile
                region_profile = profile[region_start_idx:region_end_idx + 1]
                # Score region
                score = _score_region(region_profile)
                # Initialize region
                region = Region(
                    chrom,
                    FeatureLocation(region_start, region_end),
                    id="{}{}".format(label, counter + 1),
                    type=type,
                    score=score,
                    qualifiers={
                        "source": source
                    },
                    profile=region_profile
                )
                # Add region to regions
                regions.append(region)
                # Reset indices and increase counter
                region_start_idx = None
                region_end_idx = None
                counter += 1

    # Started a region but didn't finish it (scores
    # remained above the threshold until the end)
    if region_start_idx is not None:
        # Initialize region coordinates
        region_start = start + region_start_idx
        if region_end_idx is None:
            region_end_idx = len(profile) - 1
        region_end = start + region_end_idx
        # Get region profile
        region_profile = profile[region_start_idx:region_end_idx + 1]
        # Score region
        score = _score_region(region_profile)
        # Initialize region
        region = Region(
            chrom,
            FeatureLocation(region_start, region_end),
            id="{}{}".format(label, counter + 1),
            type=type,
            score=score,
            qualifiers={
                "source": source
            },
            profile=region_profile
        )
        # Add region to regions
        regions.append(region)

    # Combine regions into larger regions that still
    # score above threshold
    if stitch:
        stitched_regions = _stitch_regions(regions, profile, exons, min_score,
                                           label, type, source)

    # Only keep regions larger than min. region length
    # that still score above threshold
    for region in sorted(regions + stitched_regions,
                         key=lambda x: x.end - x.start + 1, reverse=True):
        # Initialize
        overlap = False
        # Skip short regions
        if region.end - region.start + 1 < min_length:
            continue
        # For each region...
        for fine_region in fine_regions:
            # If regions overlap...
            if region.overlap(fine_region):
                overlap = True
                break
        # Skip overlapping regions
        if not overlap:
            fine_regions.append(region)

    # Sort regions by genomic coordinates
    fine_regions.sort(key=lambda x: x.start)

    return fine_regions


def _score_region(profile):
    """
    Compute and return score of a region defined by length of
    profile array.

    Code adapted from darenillas "PhastCons.py"
    Original function: "_score_region"
    URL: https://github.com/wassermanlab/GUD/blob/master/lib/ORCA/analysis/PhastCons.py
    """

    # Initialize
    score = 0

    for i in range(len(profile)):
        score += profile[i]

    # Avg. region score by length
    score /= len(profile)

    return score


def _stitch_regions(regions, profile, exons=[], min_score=0.6,
   label="Region", type="region", source=None):
   """
   Stitch regions into larger regions and return those that still
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

   # Initialize #
   stitched_regions = []

   # For each region... #
   for i in range(len(regions) - 1):
       stitch = True
       # 0-based index coords
       chrom = regions[i].chrom
       A_start = regions[i].start
       A_end = regions[i].end
       A_start_idx = regions[i].start - regions[0].start
       A_end_idx = regions[i].end - regions[i].start + A_start_idx
       A_region_id = i + 1
       # For each other region... #
       for j in range(i + 1, len(regions)):
           # 0-based index coords
           B_start = regions[j].start
           B_end = regions[j].end
           B_start_idx = regions[j].start - regions[0].start
           B_end_idx = regions[j].end - regions[j].start + B_start_idx
           B_region_id = j + 1
           # If filtering exons, don't combine regions
           # separated by coding exons. E.g.:
           #
           # Region : ----AAAA-----BBBB----
           # exon 1 : XXX------------------ => exon_between = False
           # exon 2 : ------------------XXX => exon_between = False
           # exon 3 : ---------XXX--------- => exon_between = True
           # exon 4 : ---XXX--------------- => exon_between = False
           # exon 5 : ---------------XXX--- => exon_between = False
           # exon 6 : ------XXX------------ => exon_between = True
           # exon 7 : ------------XXX------ => exon_between = True
           if exons:
               # Initialize
               exon_between = False
               # For each exon...
               for exon in exons:
                   exon_start = exon.start
                   exon_end = exon.end
                   # exons 2, 5
                   if exon_start > B_start:
                       continue
                   # exons 1, 4
                   if exon_end < A_end:
                       continue
                   # exons 3, 6, 7
                   if exon_end >= A_end and exon_start <= B_start:
                       exon_between = True
                       break
               # Don't combine regions if separated by exon
               if exon_between:
                   stitch = False
           # Combine regions
           if stitch:
               # Initialize
               region_profile = profile[A_start_idx:B_end_idx+1]
               score = _score_region(region_profile)
               # If combined region still scores above
               # threshold...
               if score >= min_score:
                   # Initialize combined region
                   stitched_region = Region(
                       chrom,
                       FeatureLocation(A_start, B_end),
                       id = "{}{}-{}".format(label, A_region_id, B_region_id),
                       type = type,
                       score = score,
                       qualifiers = {
                           "source" : source
                       },
                       profile = region_profile
                   )
                   # Add region to regions
                   stitched_regions.append(stitched_region)
               # Don't combine regions if scores below
               # threshold
               else:
                   stitch = False
           # Cannot combine any further
           if not stitch: break

   return stitched_regions
