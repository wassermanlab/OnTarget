import traceback

import requests
from flask import request

from flaskr.Endpoints.name2sample import name2sample


def gene2region(remote_host, database, genes, samples, mode):
    URL = remote_host + '/api/v1/' + database + '/genes'
    results = {}
    if genes is not None:
        try:
            results['results'] = {}
            for gene in genes:
                # Note we may have to update to accomodate multiple names
                data = requests.get(url=URL, params={'names':gene}).json()
                list = []
                while "next" in data:
                    list.extend(data['results'])
                    data = requests.get(url=data['next']).json()
                list.extend(data['results'])
                results['results'][gene] = get_gene_region(list, database, remote_host, samples=samples, limit_by=mode)
        except:
            traceback.print_exc()
            results["Error"] = "Could not find Genes"
    return results


def get_gene_region(gene, database, remote_host, samples=[], limit_by="tad"):
    """
    Delimits a region for the given gene based
    on TAD boundaries, nearby genes, or +/- N
    kb.
    """
    # Get chromosome sizes
    chromURL = remote_host + '/api/v1/' + database + '/chroms'
    try:
        chromData = requests.get(url=chromURL).json()
    except:
        print("error connecting to GUD")
    chrom_sizes = {}
    for item in chromData['results']:
        chrom_sizes[item['chrom']] = item['size']
    # Get all genes with the given name
    genes = gene
    # If genes are not valid...
    if not genes:
        raise ValueError(
            "Gene \"%s\" is not valid!!!" % gene
        )

    gene_start = min([item['start'] for item in genes])
    gene_end = max([item['end'] for item in genes])
    chrom = genes[0]['chrom']
    region_start = 0
    region_end = chrom_sizes[chrom]
    # If delimit by distance...
    if limit_by.isdigit():
        region_start = gene_start - int(limit_by) * 1000
        region_end = gene_end + int(limit_by) * 1000
    # # Instead, if delimit by TADs...
    elif limit_by == "tad":
        region_start, region_end =\
            get_region_coordinates_by_tad(
                gene,
                chrom,
                gene_start,
                gene_end,
                region_start,
                region_end,
                samples
            )
    #
    # Instead, if delimit by genes...
    elif limit_by == "gene":
        region_start, region_end = get_region_coordinates_by_gene(genes, region_start, region_end, gene_start, gene_end,
                                                                  URL=remote_host + '/api/v1/' + database + '/genes')

    return {'start': region_start, 'end': region_end}


def get_region_coordinates_by_gene(genes, region_start, region_end, gene_start, gene_end, URL):
    previous_gene = requests.get(url=URL, params={'uids': min([g["qualifiers"]["uid"] for g in genes]) - 1}).json()
    previous_genes = requests.get(url=URL, params={'names': previous_gene["qualifiers"]["name"]}).json()
    region_start = get_coord(genes, previous_genes, "start", URL)
    next_gene = requests.get(url=URL, params={'uids': max([g["qualifiers"]["uid"] for g in genes]) + 1}).json()
    next_genes = requests.get(url=URL, params={'names': next_gene["qualifiers"]["name"]}).json()
    region_end = get_coord(genes, next_genes, "end", URL)
    return region_start, region_end


def get_coord(A, B, coord, URL):
    genes_overlap = False
    for a in A:
        for b in B:
            if overlap(a, b):
                genes_overlap = True
    if genes_overlap:
        if coord == "start":
            previous_gene = requests.get(url=URL, params={'uids': min([b["qualifiers"]["uid"] for b in B]) - 1}).json()
            previous_genes = requests.get(url=URL, params={'names': previous_gene["qualifiers"]["name"]}).json()
            return get_coord(A, previous_genes, coord, URL)
        else:
            next_gene = requests.get(url=URL, params={'uids': min([b["qualifiers"]["uid"] for b in B]) + 1}).json()
            next_genes = requests.get(url=URL, params={'names': next_gene["qualifiers"]["name"]}).json()
            return get_coord(A, next_genes, coord, URL)
    if coord == "start":
        return max([b["end"] for b in B])
    else:
        return min([b["start"] for b in B])


def overlap(self, feat):
    """
    Returns if feat overlaps this one or not.
    """
    try:
        if self.chrom == feat.chrom and self.start < feat.end and self.end > feat.start:
            return True
        else:
            return False
    except:
        raise ValueError(
            "Could not calculate overlap!"
        )


def get_region_coordinates_by_tad(gene, chrom, gene_start, gene_end, region_start, region_end, remote_host, database, samples=[]):
    tads = []
    sampleIDs = []
    if samples:
        tads = get_tads(
            remote_host + '/api/v1/' + database + '/tads',
            chrom,
            gene_start,
            gene_end,
            region_start,
            region_end,
            samples
        )
    # TODO: add TSS function once done in GUD
    # TSS endpoint not functional in GUD yet
    # if not tads:
    #     for tss in TSS.select_by_gene(
    #         session,
    #         gene,
    #         as_genomic_feature=True
    #     ):
    #         sampleIDs +=\
    #             tss.qualifiers["sampleIDs"]
    #     # Get samples
    #     samples = []
    #     for obj in name2sample(remote_host, database)['results']:
    #         if obj['uid'] in sampleIDs:
    #             samples.append(obj)
    #
    #     tads = get_tads(
    #         chrom,
    #         gene_start,
    #         gene_end,
    #         region_start,
    #         region_end,
    #         [s.name for s in samples]
    #     )
    if not tads:
        tads = get_tads(
            remote_host + '/api/v1/' + database + '/tads',
            chrom,
            gene_start,
            gene_end,
            region_start,
            region_end
        )
    for t in tads:
        # Get upstream start closest
        # to the gene's start
        if gene_start >= t.start > region_start:
            region_start = t.start
        # Get downstream end closest
        # to the gene's end
        if gene_end <= t.end < region_end:
            region_end = t.end

    return region_start, region_end


def get_tads(URL, chrom, gene_start,
             gene_end, region_start, region_end,
             samples=[]):
    # Initialize
    encompassing_tads = []

    # Get TADs
    if len(samples) != 0:
        tads = requests.get(url=URL, params={'chrom': chrom, 'start':region_start,
                                         'end':region_end,'location':'overlapping','samples':','.join(map(str, samples))}).json()
    else:
        tads = requests.get(url=URL, params={'chrom': chrom, 'start':region_start,
                                         'end':region_end,'location':'within'}).json()

    for t in tads:
        if t.end >= gene_end and \
                t.start <= gene_start:
            encompassing_tads.append(t)

    if encompassing_tads:
        return encompassing_tads

    return tads
