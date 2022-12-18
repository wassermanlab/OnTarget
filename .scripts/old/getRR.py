
#
# TODO:
# 1) Handle multiple entries for the same gene
#    In the UCSC DB, a single gene may have multiple records corresponding
#    to different transcripts.
#    a) take longest transcript?
#    b) combine transcripts/exons?
#    c) treat all (protein coding) exons individually regardless of
#       overlap? - DONE
# 2) Fix comparison to other genes in the region cut down function.
#    a) compare 'name2' not 'name'? - DONE
# 3) Limit histone marks to H3K27Ac - DONE
# 4) Add cell or tissue type filter - DONE
# 5) Add expression level tables to FANTOM5 enhancers for selection by
#    cell/tissue type. This may be done for some future version where the
#    feature score is used in the overall scoring scheme
# 6) Use CDS start/end as gene boundaries for region cut down rather than
#    transcription start/end?
# 7) Used actual phastCons scores at each nucleotide rather than binsary
#    conserved regions - DONE. Now there is a choice of which to use.
# 8) Make region cutdown on neighbouring gene boundaries optional - DONE
#

import os
import sys
import warnings
import argparse
import csv

sys.path.insert(0, "lib")

from array import array

import Bio
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.motifs.jaspar.db import JASPAR5

import numpy as np
from pybedtools import BedTool

from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session

from GUD.bin_range import BinRange
from GUD.ORM.ref_gene import RefGene
#from GUD.ORM.tad_hESC_total_combined_domain_LO import TadhESCTotalCombinedDomainLO
#from GUD.ORM.tad_IMR90_total_combined_domain_LO import TadIMR90TotalCombinedDomainLO
#from GUD.ORM.dnase_cluster import DnaseCluster
#from GUD.ORM.manta_tfbs import MantaTfbs
#from GUD.ORM.tfbs_cluster import TfbsCluster
#from GUD.ORM.fantom5 import Enhancer
#from GUD.ORM.vista_enhancer import VistaEnhancer
#from GUD.ORM.decres import Decres
from GUD.ORM.conservation import Conservation
from GUD.ORM.dna_accessibility import DnaAccessibility
from GUD.ORM.histone_modification import HistoneModification
from GUD.ORM.repeat_mask import RepeatMask
from GUD.ORM.tad import Tad
from GUD.ORM.tf_binding import TFBinding
#from GUD.ORM.gro_seq import GROseq


####################
# ORIOL's comments:
# Remove ORCAtk, use only GUD!
####################
#from ORCA.analysis.PhastCons import PhastCons

#from connection import Database
#from pybedtools import BedTool
#import MySQLdb

####################
# ORIOL's comments:
# All these encoded parameters should be in a "config.ini" file and imported
# using ConfigParser; or alternatively, in an __init__.py file.
#
# Example:
# # Get module path #
# module_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..")
# # Append module to path #
# sys.path.append(module_path)
# # Read configuration file #
# config = ConfigParser.ConfigParser()
# config_file = os.path.join(module_path, "config.ini")
# config.read(config_file)
####################
# Defaults/ constants
GUD_DB              = 'mm10_old'
GUD_USER            = 'gud_r'
GUD_HOST            = 'gud.cmmt.ubc.ca'
HG19_FASTA_PATH     = '/home/oriol/OnTarget/data/genomes/mm10'
#CONSERVATION_TRACK  = 'phastCons100way'
MIN_CONSERVATION    = 0.6 # i.e. 60%
MIN_CONS_LENGTH     = 10
NO_TAD_FLANK_SIZE   = 100000    # If there are no TADs that can be used to
                                # defined the gene enhancer search region,
                                # then use this amount of flanking sequence.
#ENHANCER_RELATIVE_THRESHOLD = 0.7  # use relative threshold in range 0 - 1
#ENHANCER_ABSOLUTE_THRESHOLD = 3.5  # use absolute threshold based on some
                                    # (somewhat arbitrary) fraction of the
                                    # maximum sum of the weights
TOP_SCORE_PERCENTILE = 0.05         # use a top score percentile to determine
                                    # the score threshold used to call
                                    # enhancers
JASPAR_TFBS_SCORE_THRESHOLD = '85%'
HISTONE_TYPES = ['H3K4me1', 'H3K4me3', 'H3K9ac', 'H3K27ac', 'H3K36me3']
# XXX should we also include 'LINE?', 'SINE?' and 'LTR?'
REPEAT_CLASSES = ['LINE', 'SINE', 'LTR']
MIN_ENHANCER_LENGTH = 100

EVIDENCE_TYPE_DISCRETE      = 1
EVIDENCE_TYPE_CONTINUOUS    = 2

EVIDENCE_TYPES = {
    'conservation_profile'  : EVIDENCE_TYPE_CONTINUOUS,
    'conserved_regions'     : EVIDENCE_TYPE_DISCRETE,
    #'decres_enhancers'      : EVIDENCE_TYPE_DISCRETE,
    'dna_accessibility'     : EVIDENCE_TYPE_DISCRETE,
    #'fantom5_enhancers'     : EVIDENCE_TYPE_DISCRETE,
    'histone_modification'  : EVIDENCE_TYPE_DISCRETE,
    'tf_binding'            : EVIDENCE_TYPE_DISCRETE,
    #'manta_tfbs'           : EVIDENCE_TYPE_DISCRETE,
    #'chrom_hmmm'           : EVIDENCE_TYPE_DISCRETE,
    #'gro_seq'               : EVIDENCE_TYPE_DISCRETE
}

# Evidence weights. For discrete feature types this is the weight given to
# the entire region computed for a feature of this type. For continuous
# (per nucleotide) feature types, this is the maximum possible score that an
# individual nucleotide can be assigned but NOTE this ASSUMES a scale of
# 0 to max. score(!!!) as the weight is used to compute the enhancer's overall
# score from 0 to 1.
EVIDENCE_WEIGHTS = {
    'conserved_regions'     : 1.0,
    'conservation_profile'  : 1.0,
    #'decres_enhancers'      : 0.3,
    'dna_accessibility'     : 1.0,
#    'fantom5_enhancers'     : 1,
    'histone_modification'  : 1.0,
    'tf_binding'            : 1.0,
#    'gro_seq'               : 1
    #'manta_tfbs'           : 0.5,
    #'chrom_hmmm'           : 0.2,
}


####################
# ORIOL's comments:
# CHROMS and CHROM_LENGTHS (not CHROM_LENGTHS_HG19) need to be fetch from MySQL
####################
# cat /space/data/fasta/mm10/mm10.chrom.sizes | perl -e '@chrom_names=();while(<>){chomp;($chrom,$size)=split("\t",$_);push(@chrom_names,"\"$chrom\"");}print join(", ",@chrom_names),"\n";'
CHROMS = ["chr1", "chr2", "chrX", "chr3", "chr4", "chr5", "chr6", "chr7", "chr10", "chr8", "chr14", "chr9", "chr11", "chr13", "chr12", "chr15", "chr16", "chr17", "chrY", "chr18", "chr19", "chrM"]
# cat /space/data/fasta/mm10/mm10.chrom.sizes | perl -e 'while(<>){chomp;($chrom,$size)=split("\t",$_);print("\"$chrom\" : $size,\n");}'
CHROM_LENGTHS_HG19 = {
    "chr1" : 195471971,
    "chr2" : 182113224,
    "chrX" : 171031299,
    "chr3" : 160039680,
    "chr4" : 156508116,
    "chr5" : 151834684,
    "chr6" : 149736546,
    "chr7" : 145441459,
    "chr10" : 130694993,
    "chr8" : 129401213,
    "chr14" : 124902244,
    "chr9" : 124595110,
    "chr11" : 122082543,
    "chr13" : 120421639,
    "chr12" : 120129022,
    "chr15" : 104043685,
    "chr16" : 98207768,
    "chr17" : 94987271,
    "chrY" : 91744698,
    "chr18" : 90702639,
    "chr19" : 61431566,
    "chrM" : 16299
}

class Region(SeqFeature):
    """Implements a Region object.
    
    Code adapted from darenillas "conserved_region.py"
    Original class: "ConservedRegion"
    URL: https://github.com/wassermanlab/GUD/blob/master/lib/ORCA/conserved_region.py
    
    """

    def __init__(self, location=None, type='region', id='<unknown id>',
                 qualifiers=None, score=None, profile=None, ref_seq=None):

        if score is not None:
            self.score = float(score)

        if (profile is not None
            and not isinstance(profile, array)):

            raise ValueError(
                "Profile not None and is not an array"
            )

        self.profile = profile

        super(Region, self).__init__(location, type=type, id=id, qualifiers=qualifiers)

    def __str__(self):
        return "{}\t{}\t{}\t{:.3f}".format(self.id, self.start, self.end, self.score)

    def __repr__(self):
        return "{}({}, id={}, score={})".format(self.type, self.location, self.id, self.score)

    @property
    def start(self):
        return self.location.start

    @property
    def end(self):
        return self.location.end

    @property
    def start_1_based(self):
        return self.start + 1

    @property
    def end_1_based(self):
        return self.end

    def overlap(self, feat):
        """If feature overlaps this region return True;
        otherwise return False

        """

        if feat:
            # Note that, for 0-based feature, if
            # feat1.end == feat2.start, feats don't
            # overlap but are book-ended. E.g.
            #
            # ---------1111111---------
            # -------22222------------- => True
            # ----------22222---------- => True
            # -------------22222------- => True
            # -----22222--------------- => False
            # ---------------22222----- => False
            #
            if (self.start < feat.end and self.end > feat.start):
                return True

        return False

def compute_profile(start, end, features):
    """Return an array of scores from features in the given
    range.

    Code adapted from darenillas "PhastCons.py"
    Original function: "compute_conservation_profile"
    URL: https://github.com/wassermanlab/GUD/blob/master/lib/ORCA/analysis/PhastCons.py

    """

    #
    # Extract per-position scores.
    #
    scores = {}
                
    for feat in features:
        for i in range(feat.start, feat.end):
            if hasattr(feat, 'score'):
                scores.setdefault(i, feat.score)
                
            # Assign position a score of 1.0
            else:
                scores.setdefault(i, 1.0)

    #
    # Generate a profile as a continuous array data.
    # Fill gaps with 0 values.
    #        
    profile = array('f')
            
    for i in range(start, end):
        if i in scores:
            profile.append(scores[i])
        
        else:
            profile.append(0.0)
        
    return profile

def compute_regions(start, end, profile, exons=[],
    min_score=0.6, min_length=10, stitch=False,
    label="Region", type="region", source=None):
    """Compute regions from profile of at least min. length
    that score above min.

    """

    counter = 0
    start_idx = 0
    end_idx = end - start
    region_start_idx = None
    region_end_idx = None
    regions = []
    stitched_regions = []
#    # This should call 6 conserved regions, and merge 2-3-4, and 5-6
#    profile = array('f', [6.80556863e-01, 9.43449691e-01, 8.37105628e-01, 6.53790967e-01,
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
#                            9.14818075e-01, 8.34835039e-01, 6.40025535e-01])
#    end_idx = len(profile) - 1
#    start = start_idx + 1
#    end = end_idx + 1
#    exons = []

    for exon in exons:

        mask_start_idx = exon.start - start
        mask_end_idx = exon.end - start

        if mask_end_idx <= start_idx or mask_start_idx >= end_idx:
            continue

        if mask_start_idx < start_idx:
            mask_start_idx = start_idx

        if mask_end_idx > end_idx:
            mask_end_idx = end_idx

        profile[mask_start_idx:mask_end_idx] = array('f',
            [0.0] * (mask_end_idx - mask_start_idx))

    for i in range(len(profile)):
        if profile[i] >= min_score:

            region_end_idx = i

            if region_start_idx is None:
                region_start_idx = region_end_idx

        else:
            if region_start_idx is not None:
                region_start = start + region_start_idx
                region_end = start + region_end_idx
                region_profile = profile[region_start_idx:region_end_idx+1]

                score = score_region(region_profile)

                region = Region(
                    FeatureLocation(region_start, region_end),
                    id = "{}{}".format(label, counter + 1),
                    type = type,
                    score = score,
                    qualifiers = {
                        'source' : source,
                        'score': score
                    },
                    profile = region_profile
                )

                regions.append(region)

                region_start_idx = None
                region_end_idx = None
                counter += 1

    if region_start_idx is not None:
        # Started a region and didn't finish it (scores remained
        # above min. conservation until end of sequence)
        region_start = start + region_start_idx

        if region_end_idx is None:
            region_end_idx = len(profile) - 1

        region_end = start + region_end_idx
        region_profile = profile[region_start_idx:region_end_idx+1]

        score = score_region(region_profile)

        region = Region(
            FeatureLocation(region_start, region_end),
            id = "{}{}".format(label, counter + 1),
            type = type,
            score = score,
            qualifiers = {
                'source' : source,
                'score': score
            },
            profile = region_profile
        )

        regions.append(region)

    if stitch:
        #
        # Combine seed regions into larger regions that still score
        # above min. conservation
        #
        stitched_regions = stitch_regions(regions, profile, exons, min_score,
            label, type, source)

    #
    # Keep only regions >= min. region length that still score
    # above min. conservation.
    #
    fine_regions = []
    
    for region in sorted(regions + stitched_regions,
        key=lambda x: x.end - x.start + 1, reverse=True):

        # Too small
        if region.end - region.start + 1 < min_length:
            continue

        overlap = False

        for fine_region in fine_regions:

            if region.overlap(fine_region):
                overlap = True
                break

        if not overlap:
            fine_regions.append(region)

    fine_regions.sort(key=lambda x: x.start)
    
    return fine_regions

def stitch_regions(regions, profile, exons=[], min_score=0.6,
    label="Region", type="region", source=None):
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

    stitched_regions = []

    for i in range(len(regions) - 1):

        stitch = True

        # 0-based index coords
        A_start = regions[i].start
        A_end = regions[i].end
        A_start_idx = regions[i].start - regions[0].start
        A_end_idx = regions[i].end - regions[i].start + A_start_idx
        A_region_id = i + 1

        for j in range(i + 1, len(regions)):

            # 0-based index coords
            B_start = regions[j].start
            B_end = regions[j].end
            B_start_idx = regions[j].start - regions[0].start
            B_end_idx = regions[j].end - regions[j].start + B_start_idx
            B_region_id = j + 1

            if exons:
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
                #
                exon_between = False

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

                if exon_between:
                    stitch = False

            if stitch:

                region_profile = profile[A_start_idx:B_end_idx+1]
                
                score = score_region(region_profile)

                # Combine regions
                if score >= min_score:
                    stitched_region = Region(
                        FeatureLocation(A_start, B_end),
                        id = "{}{}-{}".format(label, A_region_id, B_region_id),
                        type = type,
                        score = score,
                        qualifiers = {
                            'source' : source,
                            'score': score
                        },
                        profile = region_profile
                    )

                    # Append combined region, length, and individual regions
                    stitched_regions.append(stitched_region)

                else:
                    stitch = False

            # Cannot combine any further
            if not stitch:
                break

    return stitched_regions

def score_region(profile):
        """Compute and return score of a region defined by length of
        profile array.

        Code adapted from darenillas "PhastCons.py"
        Original function: "_score_region"
        URL: https://github.com/wassermanlab/GUD/blob/master/lib/ORCA/analysis/PhastCons.py

        """

        score = 0

        for i in range(len(profile)):
            score += profile[i]

        score /= len(profile)

        return score

#####################
## New function to start by if a user gives a region, rather than a gene
#####################
#
def start_by_region(ud_region):
    
    chrom = None

    #split the coordinate string: expected chr#:int-int
    templist = ud_region.split(":",1) #get the chrom
    chrom = templist[0]
    #check if its a valid chromosome
    if chrom not in CHROM_LENGTHS_HG19:
        warnings.warn("Not a valid chromosome")
        return None

    chrom_end = CHROM_LENGTHS_HG19[chrom]

    templist = templist[1].split("-",1) #now list is made of a start and end coordinate
    gene_start = int(templist[0])
    gene_end = int(templist[1])

    #sanity checks
    if gene_start > gene_end:
        warnings.warn("Start coordinate cannot be larger than end coordinate.")
        return None
    if gene_start < 0:
        warnings.warn("Invalid start coordinate. Search space will default to 0")
        gene_start = 0
    if gene_end > chrom_end:
        warnings.warn("Search end is outside of chromosome space. Search space will default to chromosome end")
        gene_end = chrom_end

    # This is how Dave originally stored feature data, so I'm keeping it in that form for now...
    gene_region = SeqFeature(FeatureLocation(gene_start, gene_end))

    #
    # Retrieve all genes in the gene region in order to make the exon
    # mask. This is copied from the previous code as it is needed
    # downstream. 
    #
    all_genes = RefGene.select_overlapping_range(session, chrom, int(gene_region.location.start), int(gene_region.location.end))

    if verbose:
        print("\nAll genes:")
        for g in all_genes: print(g)

    #will have to return the correct sets of info: chrom, gene_region, all_genes
    return chrom, gene_region, all_genes


#########################
# Throwing the first part of Dave's original code here.
#########################
def start_by_gene(gene_symbol, do_gene_cutdown=False):

    # Fetch all genes (transcripts) with given common (HGNC Symbol) name.
    genes = RefGene.select_by_names(session, [gene_symbol])

    if verbose:
        print("\nGenes:")
        for g in genes: print(g)

    ####################
    # ORIOL's comments:
    # I don't understand this looping thing
    ####################
    # Loop through the genes
    chrom = None
    min_gene_start = 999999999
    max_gene_end   = 0
    for gene in genes:
        # ignore scaffolds etc.
        if gene.chrom in CHROMS:
            if not chrom:
                chrom = gene.chrom
            else:
                if gene.chrom != chrom:
                    sys.exit("Gene/transcript '{}' chromosome '{}' mismatches chromome from previously fetched gene/transcript chromosome '{}'\n".format(gene.name2, gene.chrom, chrom))

            if gene.start < min_gene_start:
                min_gene_start = gene.start

            if gene.end > max_gene_end:
                max_gene_end = gene.end

    if chrom == None:
        warnings.warn("Gene {} was not found on a characterized chromosome".format(gene_symbol))
        return None

    ####################
    # ORIOL's comments:
    # I modified the TAD search:
    # . find an encompassing TAD for the given cell/tissue
    # . if failed, find an overlapping TAD for the given cell/tissue
    # . if failed, find any encompassing TAD
    # . if failed, find any overlapping TAD
    ####################

    tad_feats = []

    # if no_tads:
    #     tad_feats = []
    
    # else:
    #     tad_feats = Tad.select_encompasing_range(session, chrom, min_gene_start,
    #         max_gene_end, tissue=cell_or_tissue)
    #     # If no TADs encompassing the gene are found, try overlapping instead
    #     if not tad_feats:
    #         tad_feats = Tad.select_overlapping_range(session, chrom, min_gene_start,
    #             max_gene_end, tissue=cell_or_tissue)

    #     if not tad_feats and cell_or_tissue:
    #         #
    #         # Could not get TADs for the given cell or tissue types. In this case
    #         # just use any available TADs.
    #         #
    #         warnings.warn("\nTAD Features for cell_or_tissue=%s could not be found!\n\tSwitching to cell_or_tissue=None\n" % cell_or_tissue)
    #         tad_feats = Tad.select_encompasing_range(session, chrom, min_gene_start,
    #             max_gene_end, tissue=None)
    #         # If no TADs encompassing the gene are found, try overlapping instead
    #         if not tad_feats: 
    #             tad_feats = Tad.select_overlapping_range(session, chrom, min_gene_start,
    #                 max_gene_end, tissue=None)

    #     if filter_by:
    #         filtered_feats = []

    #         for feat in tad_feats:
    #             if feat.cell_or_tissue.count(filter_by) == 0:
    #                 filtered_feats.append(feat)

    #     tad_feats = filtered_feats
      
    #
    # If TAD data, use TAD boundaries as gene region; otherwise, 
    # gene region = gene body +/- NO_TAD_FLANK_SIZE bp
    #
    gene_region_start = 0
    gene_region_end = CHROM_LENGTHS_HG19[chrom]

    if tad_feats:
        
        if verbose:
            print("\nTAD Features:")
                
        for feat in tad_feats:
            if verbose:
                print(feat)
    
            if feat.start <= min_gene_start and feat.start > gene_region_start:
                gene_region_start = feat.start

            if feat.end >= max_gene_end and feat.end < gene_region_end:
                gene_region_end = feat.end
    
    else:
        if (min_gene_start - NO_TAD_FLANK_SIZE) > gene_region_start:
            gene_region_start = min_gene_start - NO_TAD_FLANK_SIZE
        if (max_gene_end + NO_TAD_FLANK_SIZE) < gene_region_end:
            gene_region_end = max_gene_end + NO_TAD_FLANK_SIZE
    
    if verbose:
        gene_region_len = gene_region_end - gene_region_start
        print("\nGene region: {}:{}-{}\tsize={}".format(chrom, gene_region_start, gene_region_end, gene_region_len))

    gene_region = SeqFeature(FeatureLocation(gene_region_start, gene_region_end))

    #
    # Retrieve all genes in the gene region in order to cutdown the search
    # region further by the 5' end of such genes.
    # Note: this will also return the gene of interest
    #
    all_genes = RefGene.select_overlapping_range(session, chrom,
        int(gene_region.location.start), int(gene_region.location.end))

    if verbose:
        print("\nAll genes:")
        for g in all_genes: print(g)

    other_genes = filter_for_other_genes(gene_symbol, all_genes)

    if verbose:
        print("\nOther genes:")
        for g in other_genes: print(g)

    # If we specificied cutting down the region further by neighbouring genes
    # OR if the gene of interest wasn't complelely contained by TADs then
    # cutdown by neighbouring genes.
    #if (do_gene_cutdown or region_start is None or region_end is None):
    if do_gene_cutdown:
        which_side = 'both'
        #which_side = None
        #if do_gene_cutdown or (region_start is None and region_end is None):
        #    which_side = 'both'
        #elif region_start == 0:
        #    which_side = 'left'
        #elif region_end == CHROM_LENGTHS_HG19[chrom]:
        #    which_side = 'right'

        gene_region = cutdown_region_by_neighbouring_genes(
            gene_region, min_gene_start, max_gene_end, other_genes,
            which=which_side)

        if verbose:
            print("\nGene region after cutdown: {}:{}-{}\tsize={}".format(chrom, gene_region.location.start, gene_region.location.end, gene_region.location.end - gene_region.location.start))

    #will eventually return the chrom, gene_region, and all_genes
    return chrom, gene_region, all_genes 

##def compute_gene_ontarget_enhancers(
##    gene_symbol, out_bed_file, out_track_file,
##    out_tf_file=None,
##    out_scores_file=None,
##    out_unmasked_scores_file=None,
##    do_gene_cutdown=False,
##    top_score_percentile=TOP_SCORE_PERCENTILE):

#new DEF statement, replacing gene_symbol with regions and genes to avoid code duplication
def compute_RRs(chrom, gene_region, all_genes, out_bed_file,
    out_track_file, out_tf_file=None, out_scores_file=None,
    out_unmasked_scores_file=None, do_gene_cutdown=False,
    top_score_percentile=TOP_SCORE_PERCENTILE,
    min_enhancer_length=MIN_ENHANCER_LENGTH):

##    # Fetch all genes (transcripts) with given common (HGNC Symbol) name.
##    genes = RefGene.select_by_names(session, [gene_symbol])
##
##    if verbose:
##        print "\nGenes:"
##        for g in genes:
##            print g
##
##    # Loop through the genes
##    chrom = None
##    min_gene_start = 999999999
##    max_gene_end   = 0
##    for gene in genes:
##        # ignore scaffolds etc.
##        if gene.chrom in CHROMS:
##            if not chrom:
##                chrom = gene.chrom
##            else:
##                if gene.chrom != chrom:
##                    sys.exit("Gene/transcript '{}' chromosome '{}' mismatches chromome from previously fetched gene/transcript chromosome '{}'\n".format(gene.name2, gene.chrom, chrom))
##
##            if gene.start < min_gene_start:
##                min_gene_start = gene.start
##
##            if gene.end > max_gene_end:
##                max_gene_end = gene.end
##
##    if chrom == None:
##        warnings.warn("Gene {} was not found on a characterized chromosome".format(gene_symbol))
##        return None
##
##    tads = Tad.select_overlapping_range(session, chrom, min_gene_start,
##                max_gene_end, tissue=cell_or_tissue)
##
##    if not tads and cell_or_tissue:
##        #
##        # Could not get TADs for the given cell or tissue types. In this case
##        # just use any available TADs.
##        #
##        warnings.warn("No TADs defined for the given cell/tissue types; using TADs for all cell/tissue types")
##
##        tads = Tad.select_overlapping_range(session, chrom, min_gene_start,
##                    max_gene_end, tissue=None)
##
##    gene_region_start = min_gene_start - NO_TAD_FLANK_SIZE
##    gene_region_end = max_gene_end + NO_TAD_FLANK_SIZE
##
##    if gene_region_start < 0:
##        gene_region_start = 0
##
##    if gene_region_end > CHROM_LENGTHS_HG19[chrom]:
##        gene_region_end = CHROM_LENGTHS_HG19[chrom]
##
##    if tads:
##        if verbose:
##            print "\nTADs:"
##            for tad in tads:
##                print tad
##
##        for tad in tads:
##            if tad.start <= min_gene_start and tad.start > gene_region_start:
##                gene_region_start = tad.start
##
##            if tad.end >= max_gene_end and tad.end < gene_region_end:
##                gene_region_end = tad.end
##    #else:
##    #    raise ValueError("No TADs defined in gene region!")
##
##    if verbose:
##        gene_region_len = gene_region_end - gene_region_start
##        print "\nGene region: {}:{}-{}\tsize={}".format(
##            chrom, gene_region_start, gene_region_end, gene_region_len)
##
##    gene_region = SeqFeature(FeatureLocation(gene_region_start,
##                                             gene_region_end))
##
##    #
##    # Retrieve all genes in the gene region in order to cutdown the search
##    # region further by the 5' end of such genes.
##    # Note: this will also return the gene of interest
##    #
##    all_genes = RefGene.select_overlapping_range(
##                    session, chrom, gene_region.location.start,
##                    gene_region.location.end)
##
##    if verbose:
##        print "\nAll genes:"
##        for g in all_genes:
##            print g
##
##    other_genes = filter_for_other_genes(gene_symbol, all_genes)
##
##    if verbose:
##        print "\nOther genes:"
##        for g in other_genes:
##            print g
##
##    # If we specificied cutting down the region further by neighbouring genes
##    # OR if the gene of interest wasn't complelely contained by TADs then
##    # cutdown by neighbouring genes.
##    #if (do_gene_cutdown or region_start is None or region_end is None):
##    if do_gene_cutdown:
##        which_side = 'both'
##        #which_side = None
##        #if do_gene_cutdown or (region_start is None and region_end is None):
##        #    which_side = 'both'
##        #elif region_start == 0:
##        #    which_side = 'left'
##        #elif region_end == CHROM_LENGTHS_HG19[chrom]:
##        #    which_side = 'right'
##
##        gene_region = cutdown_region_by_neighbouring_genes(
##                                gene_region, min_gene_start,
##                                max_gene_end, other_genes, which=which_side)
##
##        if verbose:
##            print "\nGene region after cutdown: {}:{}-{}\tsize={}".format(
##                chrom, gene_region.location.start, gene_region.location.end,
##                gene_region.location.end - gene_region.location.start)

    #######################################################################
    ## IMO I think that the portion below should be the entry point for compute_gene_ontarget_enchangers function.
    ## I will try placing the above portion into a separate space that only needs to take in the gene name/file
    ## In this sense, we will pass in region coordinates. I do not believe that the gene name is needed in any
    ## important part of the following lines of code. I will need to probably set a flag for gene name = None.
    #######################################################################

    #
    # Pre-compute the bins to speed up database selects rather than
    # re-computing them for each feature type fetched.
    #
    bins = BinRange().allBinsInRange(gene_region.location.start,
        gene_region.location.end)

    if verbose:
        print("Bins:")
        for b in bins: print(b)

    # Coding exons
    all_coding_exons = []
    for gene in all_genes:
        all_coding_exons.extend(gene.coding_exons)

    sequence = fetch_sequence(chrom, gene_region.location.start,
        gene_region.location.end)

    #
    # DNA accessibility features
    #
    dna_accessibility_feats = DnaAccessibility.select_by_bin_range(
        session, chrom, int(gene_region.location.start),
        int(gene_region.location.end), tissue=cell_or_tissue, bins=bins)

    # if not dna_accessibility_feats:
    #     warnings.warn("\nDNA Accessibility Features for cell_or_tissue=%s could not be found!\n\tSwitching to cell_or_tissue=None\n" % cell_or_tissue)
    #     dna_accessibility_feats = DnaAccessibility.select_by_bin_range(
    #         session, chrom, gene_region.location.start,
    #         gene_region.location.end, tissue=None, bins=bins)

    # if filter_by:
    #     filtered_feats = []

    #     for feat in dna_accessibility_feats:
    #         if feat.cell_or_tissue.count(filter_by) == 0:
    #             filtered_feats.append(feat)

    #     dna_accessibility_feats = filtered_feats

    if verbose:
        print("\nDNA Accessibility features:")
        for feat in dna_accessibility_feats: print(feat)

    #
    # Histone modification features
    #
    histone_modification_feats = HistoneModification.select_by_bin_range(
        session, chrom, int(gene_region.location.start),
        int(gene_region.location.end), histone_types=HISTONE_TYPES,
        tissue=cell_or_tissue, bins=bins)
    
    # if not histone_modification_feats:
    #     warnings.warn("\nHistone Modification Features for cell_or_tissue=%s could not be found!\n\tSwitching to cell_or_tissue=None\n" % cell_or_tissue)
    #     histone_modification_feats = HistoneModification.select_by_bin_range(
    #         session, chrom, gene_region.location.start,
    #         gene_region.location.end, histone_types=HISTONE_TYPES,
    #         tissue=None, bins=bins)

    # if filter_by:
    #     filtered_feats = []

    #     for feat in histone_modification_feats:
    #         if feat.cell_or_tissue.count(filter_by) == 0:
    #             filtered_feats.append(feat)

    #     histone_modification_feats = filtered_feats

    if verbose:
        print("\nHistone Modification Features:")
        for feat in histone_modification_feats: print(feat)

    #
    # TF-binding features
    #
    tf_binding_feats = TFBinding.select_by_bin_range(
        session, chrom, int(gene_region.location.start),
        int(gene_region.location.end), tissue=cell_or_tissue, bins=bins)
    
    # if not tf_binding_feats:
    #     warnings.warn("\nTF-Binding Features for cell_or_tissue=%s could not be found!\n\tSwitching to cell_or_tissue=None\n" % cell_or_tissue)
    #     tf_binding_feats = TFBinding.select_by_bin_range(
    #         session, chrom, gene_region.location.start,
    #         gene_region.location.end, tissue=None, bins=bins)

    # if filter_by:
    #     filtered_feats = []

    #     for feat in tf_binding_feats:
    #         if feat.cell_or_tissue.count(filter_by) == 0:
    #             filtered_feats.append(feat)

    #     tf_binding_feats = filtered_feats

    if verbose:
        print("\nTF-Binding Features:")
        for feat in tf_binding_feats: print(feat)

    #
    # Conservation features
    #
    conservation_feats = Conservation.select_by_bin_range(
        session, chrom, int(gene_region.location.start),
        int(gene_region.location.end), bins=bins)

    conservation_profile = compute_profile(int(gene_region.location.start),
        int(gene_region.location.end), conservation_feats)

    conserved_regions = compute_regions(int(gene_region.location.start),
        int(gene_region.location.end), conservation_profile, exons=all_coding_exons,
        min_score=MIN_CONSERVATION, min_length=MIN_CONS_LENGTH,
        stitch=True, label="ConservedRegion", type="conserved region",
        source=conservation_feats[0].source_name)

    if verbose:
        print("\nConserved Regions:")
        for feat in conserved_regions: print(feat)

#    if compute_conserved_regions:
#        conserved_regions = phca.compute_conserved_regions(
#            min_conservation = MIN_CONSERVATION,
#            min_length       = MIN_CONS_LENGTH,
#            filter_exons     = True
#            #flank_size       = 100
#        )
#
#        if verbose:
#            print "\nConserved regions:"
#            phca.print_conserved_regions()

#    conserved_regions = Conservation.select_by_bin_range(
#                    session, chrom, gene_region.location.start,
#                    gene_region.location.end,
#                    bins=bins)
#
#    if not conserved_regions:
#        print("no conservation")
#        #conserv = Conservation.select_by_bin_range(session, chrom, gene_region.location.start, gene_region.location.end, bins=bins)
#
#    if verbose:
#        print("\nConserved sites:")
#        for ds in conserv: print(ds)

    # XXX not using these anymore. Use the TF ChIP-seq data instead.
    #manta_tfbs = MantaTfbs.select_by_bin_range(
    #                session, chrom, gene_region.location.start,
    #                gene_region.location.end,
    #                bins=bins)
#    tf_chip_seq = TFChipSeq.select_by_bin_range(
#                    session, chrom, gene_region.location.start,
#                    gene_region.location.end, cell_or_tissue=cell_or_tissue,
#                    tf_names=tf_names, bins=bins)
#    if not tf_chip_seq:
#        print "no chip"
#        tf_chip_seq = TFChipSeq.select_by_bin_range(session, chrom, gene_region.location.start, gene_region.location.end, cell_or_tissue=None, tf_names=tf_names, bins=bins)
#    if verbose:
#        print "\nTF ChIP-seq:"
#        for tfcs in tf_chip_seq:
#            print tfcs

    # XXX or this, or both?
    #clustered_tfbs = TfbsCluster.select_by_bin_range(
    #                    session, chrom, gene_region.location.start,
    #                    gene_region.location.end,
    #                    bins=bins)

    #if verbose:
        #print "\nClustered TFBS:"
        #for ct in clustered_tfbs:
        #    print ct

#    fantom5_enhancers = Enhancer.select_by_bin_range(
#                            session, chrom, gene_region.location.start,
#                            gene_region.location.end, cell_or_tissue=cell_or_tissue, bins=bins)
#    if not fantom5_enhancers:
#        fantom5_enhancers = Enhancer.select_by_bin_range(session, chrom, gene_region.location.start, gene_region.location.end, cell_or_tissue=None, bins=bins)
#
#    if verbose:
#        print "\nFANTOM5 Enhancers:"
#        for fe in fantom5_enhancers:
#            print fe

#    gro_seq = GROseq.select_by_bin_range(session, chrom, gene_region.location.start, gene_region.location.end, cell_or_tissue=cell_or_tissue, bins=bins) 
#    
#    if not gro_seq:
#        print "no gro-seq"
#        gro_seq = GROseq.select_by_bin_range(session, chrom, gene_region.location.start, gene_region.location.end, cell_or_tissue=None, bins=bins)
#    if verbose:
#        print"\nGRO-seq:"
#        for gs in gro_seq:
#            print gs

    #vista_enhancers = VistaEnhancer.select_by_bin_range(
    #                    session, chrom, gene_region.location.start,
    #                    gene_region.location.end,
    #                    bins=bins)

    #if verbose:
        #print "\nVISTA Enhancers:"
        #for ve in vista_enhancers:
        #    print ve

    #decres_enhancers = Decres.select_by_bin_range(
    #                    session, chrom, gene_region.location.start,
    #                    gene_region.location.end,
    #                    feature_type=DECRES_ENHANCER_FEATURE_TYPE, bins=bins)

    #if verbose:
    #    print "\nDECRES Enhancers:"
    #    for de in decres_enhancers:
    #        print de

    #jaspar_tfbs = phca.compute_conserved_tfbss(
    #    motifs                      = jaspar_motifs,
    #    tfbs_score_threshold        = JASPAR_TFBS_SCORE_THRESHOLD,
    #    filter_overlapping_tfbss    = True,
    #)

    #if verbose:
        #print "\nJASPAR conserved TFBS:"
        #for jt in jaspar_tfbs:
        #    print jt

    #
    # Create the feature score matrix. Each column of the matrix
    # corresponds to the genomic position of the features and each
    # feature type is represented as a row of the matrix so that a
    # cell in the matrix is the score for that feature type at that
    # position.
    #
    #
    # Need to unleash evidence!!! 

    
    # Unwind TF binding by TF name 
    unwind_feats = {}
    for feat in tf_binding_feats:
        unwind_feats.setdefault(feat.tf_name, [])
        unwind_feats[feat.tf_name].append(feat)
    
    tf_binding_feats = []
    for tf_name in unwind_feats:
        tf_binding_feats.append(unwind_feats[tf_name])
    
    # Unwind DNA accessibility features by experiment type
    unwind_feats = {}
    for feat in dna_accessibility_feats:
        unwind_feats.setdefault(feat.experiment_type, [])
        unwind_feats[feat.experiment_type].append(feat)
    
    dna_accessibility_feats = []
    for experiment_type in unwind_feats:
        dna_accessibility_feats.append(unwind_feats[experiment_type])

    # Unwind histone modifications by histone type 
    unwind_feats = {}
    for feat in histone_modification_feats:
        unwind_feats.setdefault(feat.histone_type, [])
        unwind_feats[feat.histone_type].append(feat)
    
    histone_modification_feats = []
    for histone_type in unwind_feats:
        histone_modification_feats.append(unwind_feats[histone_type])
    
    evidence = {
#        'fantom5_enhancers'     : fantom5_enhancers,
        'tf_binding'            : tf_binding_feats,
        'dna_accessibility'     : dna_accessibility_feats,
        'histone_modification'  : histone_modification_feats,
        'conservation_profile'  : [conservation_profile]
#        'gro_seq'               : gro_seq
        #'decres_enhancers'      : decres_enhancers
        #dnase_clusters,
        #manta_tfbs,
        #clustered_tfbs,
        #jaspar_tfbs,
        #vista_enhancers,
    }

    if conserved_regions:
        evidence['conserved_regions'] = [conserved_regions]
#    else:
#        evidence[

    feature_score_matrix = create_weighted_feature_matrix(
        gene_region, evidence, EVIDENCE_WEIGHTS)

    if verbose:
        print("\nFeature scoring matrix:")
        print(feature_score_matrix)

    # The maximum possible score of any enhancer is the sum of possible
    # evidence weights
    max_evidence_score = 0
    for e in evidence.keys():
        max_evidence_score += EVIDENCE_WEIGHTS[e]
        if verbose:
            print("\nMax. evidence score: {}".format(max_evidence_score))

    #
    # Sum each column in the feature matrix to get a vector of the sum
    # of the individual feature scores at each position.
    #
    score_vector = np.sum(feature_score_matrix, axis=0)
    print(score_vector)

    if out_unmasked_scores_file:
        with open(out_unmasked_scores_file, 'w') as of:
            for score in score_vector:
                of.write("{}\n".format(score))

    if verbose:
        print("\nFeature score vector:")
        print(score_vector)

    # Create an exon mask to mask out protein coding exons
    exon_mask_vector = create_feature_mask(gene_region, all_coding_exons)

    if verbose:
        print("\nExon mask vector:")
        print(exon_mask_vector)

    repeats = RepeatMask.select_by_bin_range(
                session, chrom, int(gene_region.location.start),
                int(gene_region.location.end), repeat_classes=REPEAT_CLASSES,
                bins=bins)

    if verbose:
        print("\nRepeats:")
        for r in repeats:
            print(r)

    repeat_mask_vector = create_feature_mask(gene_region, repeats)

    if verbose:
        print("\nRepeat mask vector:")
        print(repeat_mask_vector)

    # Mask out protein coding exons with zeros
    #score_vector = np.multiply(score_vector, exon_mask_vector)
    #score_vector = np.multiply(score_vector, repeat_mask_vector)
    score_vector *=  exon_mask_vector
    score_vector *=  repeat_mask_vector

    if verbose:
        print("\nMasked feature score vector:")
        print(score_vector)

    if out_scores_file:
        with open(out_scores_file, 'w') as of:
            for score in score_vector:
                of.write("{}\n".format(score))

    #
    # Divide each value in the masked feature vector by the number of
    # features used to compute the scores.
    #
    #(height, width) = np.shape(feature_score_matrix)
    #weighted_score_vector = np.divide(score_vector.astype('float32'), height)
    #

    #
    # If we have some continuous evidence, scale each position score based on
    # the maximum possible score.
    #
    # XXX
    # Now just using conserved regions vs. conservation profile. Will need
    # improvement if we have more possibilities for continuous vs. discrete
    # evidence # types in the future!
    # XXX
    #
    if not conserved_regions:
        print("here")
        score_vector /= max_evidence_score

    # Use the top score percentile to determine the scoring threshold used to
    # define enhancer boundaries.
    sorted_vector = np.sort(score_vector)

    top_pct_pos = int(sorted_vector.size - sorted_vector.size
                      * top_score_percentile)

    threshold = sorted_vector[top_pct_pos]

#    enhancers = compute_ontarget_enhancers(
#                    gene_region, score_vector,
#                    #ENHANCER_ABSOLUTE_THRESHOLD
#                    threshold, min_enhancer_length)

    regulatory_regions = compute_regions(gene_region.location.start,
    gene_region.location.end, array('f', score_vector),
    exons=all_coding_exons, min_score=threshold, min_length=min_enhancer_length,
    stitch=False, label="RegulatoryRegion", type="enhancer",
    source="OnTarget")

    # Identify promoters
    for rr in regulatory_regions:

        for gene in all_genes:
            if rr.end < gene.txStart:
                continue

            if rr.start > gene.txStart:
                continue

            rr.type = "promoter"
            break

    for rr in regulatory_regions:
        #
        # If we have no continuous evidence, scale each enhancer score based
        # on the maximum possible score.
        # XXX
        # Now just using conserved regions vs. conservation profile. Will need
        # improvement if we have more possibilities for continuous vs.
        # discrete evidence # types in the future!
        # XXX
        #
        #
        rr.qualifiers['score'] /= max_evidence_score

        if verbose:
            #print en
            print("{} ({}) {}\t{}\tlen={}\tscore={:.3f}".format(
                rr.id, rr.type, rr.location.start, rr.location.end,
                rr.location.end - rr.location.start, rr.qualifiers['score']))

    if out_bed_file:
        write_regulatory_regions_bed_file(out_bed_file, chrom, regulatory_regions)

    if out_track_file:
        write_ucsc_track_file(
            out_track_file, chrom, gene_region,
            enhancers,
            tads=tads,
            conserved_regions=conserved_regions,
            dna_accessibility=dnase_sites,
            histone_modification=histone_modifications,
            fantom5_enhancers=fantom5_enhancers,
            #decres_enhancers=decres_enhancers,
            tf_chip_seq=tf_chip_seq
            #clustered_tfbs=clustered_tfbs,
            #manta_tfbs=manta_tfbs,
            #dnase_clusters=dnase_clusters,
            #vista_enhancers=vista_enhancers
        )

    if out_tf_file and tf_names:
        write_enhancer_tf_file(out_tf_file, chrom, enhancers, tf_names)

    return (chrom, regulatory_regions)

def filter_for_other_genes(gene_symbol, genes):
    """ Filter list of genes for genes whose symbol does not match the given
    gene symbol.
    """

    other_genes = []
    for g in genes:
        if g.name2 != gene_symbol:
            other_genes.append(g)

    return other_genes

def cutdown_region_by_neighbouring_genes(
        region, min_start, max_end, other_genes, which=None
    ):
    """ Given a larger region encompassing a smaller region (i.e. a gene
    within the larger region), cut down the larger region based on the
    proximal genes. The min_start and max_end refer to the smaller inner
    region (i.e. the gene we are encompassing). Don't cut down within this
    start/end.
    """

    #
    # NOTE: FeatureLocations are immutable so can't set start and end of
    # the passed in region directly, thus we compute the modified start
    # and end and return a new SeqFeature object.
    #
    reg_start = region.location.start
    reg_end = region.location.end

    for og in other_genes:
        #
        # Note the logic is the same regardless of whether the gene of
        # interest is a + or - strand gene.
        #
        if og.strand == '+' or og.strand == 1:
            if (og.start > reg_start and og.start < min_start):
                reg_start = og.start
            elif (og.start < reg_end and og.start > max_end):
                reg_end = og.start
        elif og.strand == '-' or og.strand == -1:
            if (og.end > reg_start and og.end < min_start):
                reg_start = og.end
            elif (og.end < reg_end and og.end > max_end):
                reg_end = og.end

    if which == 'both':
        cutdown_region = SeqFeature(FeatureLocation(reg_start, reg_end))
    elif which == 'left':
        cutdown_region = SeqFeature(FeatureLocation(reg_start,
                                                    region.location.end))
    elif which == 'right':
        cutdown_region = SeqFeature(FeatureLocation(region.location.start,
                                                    reg_end))
    else:
        cutdown_region = region

    return cutdown_region


def create_weighted_feature_matrix(region, feature_types, weights):
    """ Create a matrix of 1's and 0's where each row represents some kind
    of genomic feature type. For each position comprising an individual
    feature within a feature type, set the evidence at that position to
    1 multiplied by the corresponding weight. Set the background to 0.

    The weigths and features are dictionaries where the keys are the
    corresponding feature names. NOTE: the keys must agree between the
    weights and features dictionaries!
    """
    
    first_feat = True
    region_start = region.location.start
    region_end   = region.location.end

    # For each type of feature (e.g. DNA accessibility)
    for feat_type in feature_types:
        if feat_type:

            weight = weights[feat_type]
            
            # This iterates through every individual features
            # of that type (e.g. DNase-seq or ATAC-seq)
            for feats in feature_types[feat_type]:

                # Initialize feature vector to 0
                feat_vect = np.zeros(region_end - region_start, np.float)

                if EVIDENCE_TYPES[feat_type] == EVIDENCE_TYPE_DISCRETE:
                    print(feat_type)
                    #
                    # Discrete (interval) evidence type (region with start and end
                    # where every position in the range is assigned the weight for
                    # this evidence type.
                    #
                    for feat in feats:
                        # Set feature positions within vector to the weight of
                        # the given feature.
                        feat_vect[feat.start - region_start:
                            feat.end - region_start] = weight

                    if first_feat:
                        # Add first weighted feature vector to matrix
                        feature_matrix = feat_vect
                        first_feat = False
                    else:
                        # Add additional weighted feature vector to matrix
                        feature_matrix = np.vstack((feature_matrix, feat_vect))

                elif EVIDENCE_TYPES[feat_type] == EVIDENCE_TYPE_CONTINUOUS:
                    #
                    # Continuous evidence type. Each position has a score. This is
                    # already an array of values.
                    #
                    if first_feat:
                        # Add first feature scoring array to matrix
                        feature_matrix = feats
                        first_feat = False
                    else:
                        # Add additional weighted feature vector to matrix
                        feature_matrix = np.vstack((feature_matrix, feats))
        else:
            sys.exit("Unknown evidence type {}\n".format(
                EVIDENCE_TYPES[feat_type]))

    return feature_matrix

def create_feature_mask(region, features):
    
    region_start = region.location.start
    region_end   = region.location.end

    # Initialize exon mask vector to 1's
    region_len = region_end - region_start
    mask = np.ones(region_len, np.int)

    for f in features:
        # Check feature actually overlaps the region. It should due to the
        # way features are selected in the main program, but this is a sanity
        # check.
        if f.end > region_start and f.start < region_end:
            mask_start = f.start - region_start
            mask_end   = f.end - region_start

            # If feature only partially overlaps the region, set mask
            # coords accordingly.
            if mask_start < 0:
                mask_start = 0
            if mask_end > region_len:
                mask_end = region_len

            # Set feature position within vector to 0
            mask[mask_start : mask_end] = 0

    return mask


def compute_ontarget_enhancers(region, feat_array, threshold, min_length):
    
    region_start = region.location.start
    region_end   = region.location.end
    #print min_length

    start = 0
    end = 0
    pos = 0
    last_val = 0
    regions = []
    for val in feat_array:
        if val >= threshold and last_val < threshold:
                # Start a new region
                start = pos
                #end = pos
        elif val < threshold and last_val >= threshold:
                # End current region
                end = pos

                length = end - start
                if length >= min_length:
                    regions.append(
                        SeqFeature(
                            FeatureLocation(
                                start + region_start, end + region_start),
                            type="Enhancer",
                            qualifiers={'score' : float(sum(
                                feat_array[start:end])) / length}))

        last_val = val
        pos += 1

    if val >= threshold:
        # Ended in a region
        end = pos
        length = end - start
        if length >= min_length:
            regions.append(SeqFeature(
                FeatureLocation(start + region_start, end + region_start),
                    type="Enhancer",
                    qualifiers={'score' : float(sum(
                        feat_array[start:end])) / length}))
                
    return regions


#
# A version of this is also defined defineRR.py. We should probably put
# this in some sort of utility function module. This one differs in that
# it is a standalone function, not a method and returns a Bio.Seq.Seq object,
# not just a string sequence.
#
def fetch_sequence(chrom=None, start=None, end=None):
    '''Retrieve a sequence based on a bedfile.

    Start and end are assumed to be specified as 0-indexed half-open
    coordinates (BED standard).
    '''

    if chrom is None:
        raise ValueError("No chromosome provided")

    if start is None:
        raise ValueError("No start position provided")

    if end is None:
        raise ValueError("No end position provided")

    # Always fetch +ve strand sequence.
    bedfile = BedTool("{:s} {:d} {:d} base 0 {:s}".format(
                        chrom, start, end, '+'),
                        from_string=True)

    fasta_file = os.path.join(HG19_FASTA_PATH, "mm10.fa")

    seq = bedfile.sequence(fi=fasta_file, s=True)

    fasta = open(seq.seqfn).read()

    #to just get coordinates, we should mask the first line of the string
    seq = Seq((fasta.split("\n",2))[1], unambiguous_dna)

    return seq


def write_regulatory_regions_bed_file(out_file, chrom, regulatory_regions):

    with open(out_file, 'w') as of:
        rr_num = 1
        for rr in regulatory_regions:
            of.write("{}\t{}\t{}\tRegulatoryRegion{}\t{:.3f}\t{}\n".format(
                chrom,
                rr.location.start,
                rr.location.end,
                rr_num,
                rr.qualifiers['score'],
                rr.type))

            rr_num += 1

def write_ucsc_track_file(out_file, chrom, region, enhancers,
                tads=None,
                conserved_regions=None,
                dna_accessibility=None,
                histone_modification=None,
                fantom5_enhancers=None,
                gro_seq=None,
                decres_enhancers=None,
                tf_chip_seq=None,
                clustered_tfbs=None,
                manta_tfbs=None,
                dnase_clusters=None,
                vista_enhancers=None
            ):

    with open(out_file, 'w') as of:
        of.write("browser position {}:{}-{}\n\n".format(
            chrom, region.location.start + 1, region.location.end))

        if tads:
            of.write("\ntrack name=\"TADs\" description=\"TADs\" visibility=pack\n")
            for feat in tads:
                of.write("{}\t{}\t{}\t{}\n".format(chrom, feat.start, feat.end, feat.tissue.replace(" ", "_")))

        if conserved_regions:
            of.write("\ntrack name=\"Conserved Regions\" description=\"ORCA Conserved Regions\" useScore=1 visibility=dense\n")
            for feat in conserved_regions:
                of.write("{}\t{}\t{}\t{}\t{:d}\n".format(
                    chrom, feat.start, feat.end, feat.id,
                    int(feat.score * 1000)))

        if dna_accessibility:
            of.write("\ntrack name=\"DNA Accessibility\" description=\"DNA Accessibility\" visibility=dense\n")
            for feat in dna_accessibility:
                of.write("{}\t{}\t{}\t{}\n".format(
                    feat.chrom, feat.start, feat.end,
                    feat.cell_or_tissue.replace(" ", "_")))

        if histone_modifications:
            of.write("\ntrack name=\"Histone Modifications\" description=\"Histone Modifications\" visibility=dense\n")
            for feat in histone_modifications:
                of.write("{}\t{}\t{}\t{}\n".format(
                    feat.chrom, feat.start, feat.end,
                    feat.cell_or_tissue.replace(" ", "_")))

        if manta_tfbs:
            of.write("\ntrack name=\"MANTA TFBS\" description=\"MANTA TFBS\" useScore=1 visibility=dense\n")
            for feat in manta_tfbs:
                of.write("{}\t{}\t{}\t{}\t{:d}\t{}\n".format(
                    feat.chrom, feat.start, feat.end, feat.name,
                int(feat.score * 1000), feat.strand))

        if tf_chip_seq:
            of.write("\ntrack name=\"TF ChIP-seq\" description=\"TF ChIP-seq\" visibility=dense\n")
            for feat in tf_chip_seq:
                of.write("{}\t{}\t{}\t{}\n".format(
                    feat.chrom, feat.start, feat.end,
                    feat.tf_name.replace(" ", "_")))

        if decres_enhancers:
            of.write("\ntrack name=\"DECRES Enhancers\" description=\"DECRES Enhancers\" useScore=1 visibility=dense\n")
            for feat in decres_enhancers:
                of.write("{0}\t{1}\t{2}\t{0}:{1}-{2}\t{3:d}\t{4}\n".format(
                    feat.chrom, feat.start, feat.end, int(feat.score * 1000),
                    feat.strand))

        if fantom5_enhancers:
            # Note: score is just a count so not useful to display. Would have
            # to link to expression table for a particular sample to get actual
            # meaningful scores.
            of.write("\ntrack name=\"FANTOM5 Enhancers\" description=\"FANTOM5 Enhancers\" visibility=dense\n")
            for feat in fantom5_enhancers:
                of.write("{}\t{}\t{}\t{}\n".format(
                    feat.chrom, feat.start, feat.end))

        if enhancers:
            of.write("\ntrack name=\"OnTarget Enhancers\" description=\"OnTarget Computed Enhancers\" useScore=1 visibility=dense\n")
            for feat in enhancers:
                of.write("{0}\t{1}\t{2}\t{0}:{1}-{2}\t{3:d}\n".format(
                    chrom, feat.location.start, feat.location.end,
                    int(feat.qualifiers['score'] * 1000)))

        if dnase_clusters:
            of.write("\ntrack name=\"DNase Clusters\" description=\"DNase Clusters\" useScore=1 visibility=dense\n")
            for feat in dnase_clusters:
                of.write("{}\t{}\t{}\t{}\t{:d}\n".format(
                    feat.chrom, feat.start, feat.end, feat.name, feat.score))

        if vista_enhancers:
            of.write("\ntrack name=\"VISTA Enhancers\" description=\"VISTA Enhancers\" useScore=1 visibility=dense\n")
            for feat in vista_enhancers:
                of.write("{}\t{}\t{}\t{}\t{:d}\n".format(
                    feat.chrom, feat.start, feat.end, feat.name, feat.score))


def write_enhancer_tf_file(out_file, chrom, enhancers, tf_names):
    with open(out_file, 'w') as of:
        enhancer_num = 1
        for en in enhancers:
            tf_dict = {}
            for tf_name in tf_names:
                tf_dict[tf_name] = 0

            enhancer_tfs = TFChipSeq.select_by_bin_range(
                            session, chrom, en.location.start, en.location.end,
                            compute_bins=True)

            for etf in enhancer_tfs:
                tf_dict[etf.tf_name] = 1

            of.write("{}\tEnhancer{}\t{}\t{}\t{}\t".format(
                gene_symbol, enhancer_num, chrom, en.location.start,
                en.location.end)
            )

            for tf_name in tf_names:
                of.write("{}".format(tf_dict[tf_name]))

            of.write("\n")

            enhancer_num += 1


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Compute putative enhancer regions."
    )

    parser.add_argument(
        '-gs', '--gene-symbol',
        help="Optional: the HGNC symbol (name) of a gene."
    )

    parser.add_argument(
        '-gf', '--gene-file',
        help="Optional: file containing gene symbols (names)"
    )

    parser.add_argument(
        '-p', '--percentile', default=TOP_SCORE_PERCENTILE,
        help="Optional: Use this top percentile of scores to determine the threshold for putative enhancers."
    )

    parser.add_argument(
        '-l', '--min_length', default=MIN_ENHANCER_LENGTH,
        help="Optional: Use this as the minimum length of predicted enhancers to report."
    )

    parser.add_argument(
        '-c', '--cell-or-tissue', nargs='*', help="Optional. Limit features to those associated to the given cell lines, tissue types or other biological conditions."
    )

    parser.add_argument(
        '-f', '--filter-by', help="Optional. Filters out features associated to the given cell lines, tissue types or other biological conditions."
    )
    
    parser.add_argument(
        '-ccr', '--compute_conserved-regions', action='store_true',
        help="Optional. If specified, compute (discrete) conserved regions for weighting rather than using the (continuous) basewise conservation scores."
    )

    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help="Optional. If specified, print out detailed information about various features and intermediate calculations"
    )

    parser.add_argument(
        '-gc', '--gene-cutdown', action='store_true', help="Optional. If specified cutdown region at neighbouring genes 5' boundaries"
    )

    parser.add_argument(
        '-d', '--dir', help="Optional. If specified, all output files are written to this directory"
    )

    parser.add_argument(
        '-o', '--out-file', help="Optional. Output predicted enhancers to the specified BED formatted file."
    )

    parser.add_argument(
        '-os', '--out-scores-file', help="Optional. Output the enhancer prediction per nucleotide scores (after masking) to the given file"
    )

    parser.add_argument(
        '-ous', '--out-unmasked_scores-file', help="Optional. Output the enhancer prediction per nucleotide (unmasked) scores to the given file"
    )

    parser.add_argument(
        '-ot', '--out-track-file', help="Optional. Output UCSC browser track file"
    )

    parser.add_argument(
        '-tf', '--tf-file', help="Optional. Input list of TF names in the TF ChIP-seq dataset."
    )

    parser.add_argument(
        # Outputs a file with the TF binding for each predicted OnTarget
        # enhancer. The format of the file is:
        #   gene_symbol, enhancer#, chrom, start, end, tf_presence_bit_string
        '-otf', '--out-tf-file',
        help="Optional: optional output TF binding file"
    )
    
    ##AS OF 2018-03-02 add in flag for region start/end
    parser.add_argument(
        '-rg', '--region', help="Skips the gene step and searches for enhancers within the region. Takes in as argument form chr#;####-####"
        )
    
    #
    # # Filters:
    # # If true, then they apply.
    # #
    # parser.add_argument(
    #     '-nt', '--no-tads', action='store_true', help="Optional. Skips TAD boundaries for selecting regulatory regions."
    # )
    

    # GET ALL ARGUMENTS FROM THE SHELL
    # Does this do it by the flag name??? --> I think so but double check. And then I can remove this note.
    args = parser.parse_args()

    gene_symbol = args.gene_symbol
    gene_file = args.gene_file

    # Should this allow multiple values? YES!!!!
    cell_or_tissue = args.cell_or_tissue

    filter_by = args.filter_by

    min_enhancer_length = int(args.min_length)

    top_score_percentile = float(args.percentile)

    compute_conserved_regions = args.compute_conserved_regions

    verbose = args.verbose

    do_gene_cutdown = args.gene_cutdown

    tf_file = args.tf_file

    user_region = args.region

    out_dir = args.dir
    out_bed_file = args.out_file
    out_scores_file = args.out_scores_file
    out_unmasked_scores_file = args.out_unmasked_scores_file
    out_track_file = args.out_track_file
    out_tf_file = args.out_tf_file
    
    #
    # Boolean filters:
    #
    # no_tads = args.no_tads
    
    #############
    ## Need to change this to reflect it can also take in a region, not only a gene/gene file. TBD
    ## 2018-03-06 : included the changes
    #############

    if not gene_symbol and not gene_file and not user_region:
        sys.exit("No gene symbol or gene symbols file or region specified!")
        
    if gene_symbol and gene_file:
        sys.exit("Please specify either a gene symbol or gene file - not both!")

    if (gene_symbol or gene_file) and user_region:
        sys.exit("Please specify a gene name or a region - not both!")

    ####################
    ## I think this needs to be changed to ONTARGET, rather than GUD. I will play with this later. Need to change global vars above?
    ####################
    
    gud_connect_str ="mysql+pymysql://{}:@{}:5506/{}".format(GUD_USER, GUD_HOST, GUD_DB)
    # Set echo=True for SQLAlchemy debugging
    engine = create_engine(gud_connect_str, echo=False)
    session = Session(engine)
    
    #jdb = JASPAR5(
    #    host     = 'vm5.cmmt.ubc.ca',
    #    name     = 'JASPAR_2016',
    #    user     = 'jaspar_r',
    #    password = ''
    #)

    #jaspar_motifs = jdb.fetch_motifs(
    #    collection  = 'CORE',
    #    tax_group   = 'vertebrates',
    #    min_ic      = 8
    #)

    gene_enhancer = {}

    tf_names = []
    if tf_file:
        with open(tf_file) as tf_fh:
            for line in tf_fh.readlines():
                tf_names.append(line.rstrip())

    ############
    ## Need to add another elif for region, rather than just a gene. Can probably otherwise copy the format from the first 'if' statement
    ############

    if gene_symbol:

        chrom, gene_region, all_genes = start_by_gene(gene_symbol,
            do_gene_cutdown=do_gene_cutdown)

        if not out_bed_file:
            out_bed_file = "{}.enhancers.bed".format(gene_symbol)

        if out_dir:
            out_bed_file = os.path.join(out_dir, out_bed_file)

            if out_track_file:
                out_track_file = os.path.join(out_dir, out_track_file)
            if out_tf_file:
                out_tf_file = os.path.join(out_dir, out_tf_file)
            if out_scores_file:
                out_scores_file = os.path.join(out_dir, out_scores_file)
            if out_unmasked_scores_file:
                out_unmasked_scores_file = os.path.join(
                                    out_dir, out_unmasked_scores_file)

##        compute_gene_ontarget_enhancers(
##            gene_symbol, out_bed_file, out_track_file,
##            out_tf_file=out_tf_file,
##            out_scores_file=out_scores_file,
##            out_unmasked_scores_file=out_unmasked_scores_file,
##            do_gene_cutdown=do_gene_cutdown,
##            top_score_percentile=top_score_percentile)

        compute_RRs(chrom, gene_region, all_genes, out_bed_file, out_track_file,
                    out_tf_file=out_tf_file, out_scores_file = out_scores_file,
                    out_unmasked_scores_file = out_unmasked_scores_file,
                    top_score_percentile = top_score_percentile, min_enhancer_length = min_enhancer_length)

    elif gene_file:
        with open(gene_file, 'r') as gene_fh:
            for line in gene_fh.readlines():
                gene_symbol = line.rstrip()

                chrom, gene_region, all_genes = start_by_gene(gene_symbol,
                    do_gene_cutdown=do_gene_cutdown)

                out_bed_file = "{}.enhancers.bed".format(gene_symbol)
                out_track_file = "{}.ucsc_track.txt".format(gene_symbol)
                #out_tf_file = "{}.enhancer_tfs.txt".format(gene_symbol)
                #out_scores_file = "{}.scores.txt".format(gene_symbol)
                #out_unmaksed_scores_file = "{}.unmaksed_scores.txt".format(gene_symbol)
                if out_dir:
                    out_bed_file = os.path.join(out_dir, out_bed_file)
                    out_track_file = os.path.join(out_dir, out_track_file)
                    #out_tf_file = os.path.join(out_dir, out_tf_file)
                    #out_scores_file = os.path.join(out_dir, out_scores_file)
                    #out_unmaksed_scores_file = os.path.join(
                    #   out_dir, out_unmaksed_scores_file)

##                compute_gene_ontarget_enhancers(
##                    gene_symbol, out_bed_file, out_track_file,
##                    out_tf_file=out_tf_file,
##                    out_scores_file=out_scores_file,
##                    out_unmasked_scores_file=out_unmasked_scores_file,
##                    do_gene_cutdown=do_gene_cutdown,
##                    top_score_percentile=top_score_percentile)

                compute_RRs(chrom, gene_region, all_genes, out_bed_file, out_track_file,
                            out_tf_file=out_tf_file, out_scores_file = out_scores_file,
                            out_unmasked_scores_file = out_unmasked_scores_file,
                            top_score_percentile = top_score_percentile, min_enhancer_length = min_enhancer_length)
                
    elif user_region:
        chrom, gene_region, all_genes = start_by_region(user_region)
        #then throw these values into compute_enhancers

        #will have to check this later, it is copied from Dave's code above
        if not out_bed_file:
            out_bed_file = "{}.enhancers.bed".format(user_region)

        if out_dir:
            out_bed_file = os.path.join(out_dir, out_bed_file)

            if out_track_file:
                out_track_file = os.path.join(out_dir, out_track_file)
            if out_tf_file:
                out_tf_file = os.path.join(out_dir, out_tf_file)
            if out_scores_file:
                out_scores_file = os.path.join(out_dir, out_scores_file)
            if out_unmasked_scores_file:
                out_unmasked_scores_file = os.path.join(
                                    out_dir, out_unmasked_scores_file)
        compute_RRs(
            chrom, gene_region, all_genes, out_bed_file, out_track_file,
            out_tf_file=out_tf_file,
            out_scores_file = out_scores_file,
            out_unmasked_scores_file = out_unmasked_scores_file,
            top_score_percentile = top_score_percentile, min_enhancer_length = min_enhancer_length)
