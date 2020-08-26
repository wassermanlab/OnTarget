from array import array
from Bio.SeqFeature import FeatureLocation, SeqFeature

class GenomicFeature(SeqFeature):
    """
    Implements a Genomic Feature object based on the Biopython's Sequence
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
        start,
        end,
        score=0,
        strand=None,
        feat_type="Feature",
        feat_id="NA",
        qualifiers=None,
        profile=None
    ):

        self.chrom = chrom
        self._start = start
        self._end = end
        self.location = FeatureLocation(
            self.start,
            self.end
        )
        self.score = score
        self._strand = strand
        self.strand = self.strand_binary
        self.type = feat_type
        self.id = feat_id
        self.qualifiers = qualifiers

        if profile is not None:
            if not isinstance(profile, array):
                raise ValueError(
                    "Input profile is not an array!"
                )

        self.profile = profile

    @property
    def start(self):

        return int(self._start)

    @property
    def start_1_based(self):

        return self.start + 1

    @property
    def end(self):

        return int(self._end)

    @property
    def end_1_based(self):

        return self.end

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

    def overlaps(self, feat):
        """
        Returns if feat overlaps this one or not.
        """

        try:
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
            if self.chrom == feat.chrom and \
                    self.start < feat.end and \
                    self.end > feat.start:

                return True
        except:
            raise ValueError(
                "Could not calculate overlap!"
            )

    def __str__(self):

        return "{}\t{}\t{}\t{}\t{}\t{}".\
            format(
                self.chrom,
                self.start,  # 0-based for BED format
                self.end,
                self.id,
                self.score,
                self.strand_string
            )

    def __repr__(self):

        return "<%s(%s, %s, %s, %s, %s, %s)>" % \
            (
                self.type,
                "chrom={}".format(self.chrom),
                "start={}".format(self.start),
                "end={}".format(self.end),
                "id={}".format(self.id),
                "score={}".format(self.score),
                "strand={}".format(self.strand)
            )

    def serialize(self):
        return {
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end,
            'id': self.id,
            'score': self.score,
            'strand': self.strand,
            'qualifiers': self.qualifiers,
        }


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
               region_profile = profile[region_start_idx:region_end_idx+1]
               # Score region
               score = _score_region(region_profile)
               # Initialize region
               region = Region(
                   chrom,
                   FeatureLocation(region_start, region_end),
                   id = "{}{}".format(label, counter + 1),
                   type = type,
                   score = score,
                   qualifiers = {
                       "source" : source
                   },
                   profile = region_profile
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
       region_profile = profile[region_start_idx:region_end_idx+1]
       # Score region
       score = _score_region(region_profile)
       # Initialize region
       region = Region(
           chrom,
           FeatureLocation(region_start, region_end),
           id = "{}{}".format(label, counter + 1),
           type = type,
           score = score,
           qualifiers = {
               "source" : source
           },
           profile = region_profile
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

#-------------#
# Private     #
#-------------#

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