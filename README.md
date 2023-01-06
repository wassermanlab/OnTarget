# OnTarget
> This is a webserver to design promoters with "on target" expression patterns.

## Manifest
+ app - Webserver
+ data - Examples, genomes, liftOver chains
+ ontarget - The OnTarget module

## Requirements
OnTarget requires the following dependencies:
* [`GUD`](https://github.com/wassermanlab/GUD) with [`FuzzyWuzzy`](https://github.com/seatgeek/fuzzywuzzy), [`interval-binning`](https://interval-binning.readthedocs.io/en/latest/), [`PyMySQL`](https://pymysql.readthedocs.io/en/latest/), [`SQLAlchemy`](https://www.sqlalchemy.org/), [`SQLAlchemy-FullText-Search`](https://github.com/mengzhuo/sqlalchemy-fulltext-search), [`SQLAlchemy-Utils`](https://sqlalchemy-utils.readthedocs.io/en/latest/) 
<!-- * [`liftOver`](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) -->
* [`Python 3.x`](https://www.python.org) with:
    - [`Biopython`](https://biopython.org)
    - [`Click`](https://click.palletsprojects.com/en/8.1.x/) with [`click-option-group`](https://click-option-group.readthedocs.io/en/latest/)
    - [`flask`](https://flask.palletsprojects.com/en/1.0.x/) with [`Flask-CORS`](https://flask-cors.readthedocs.io/en/latest/) and [`Flask-RESTful`](https://flask-restful.readthedocs.io/en/latest/)
    - [`genomepy`](https://vanheeringen-lab.github.io/genomepy/)
    - [`NumPy`](https://numpy.org/)
    - [`pandas`](https://pandas.pydata.org/)
    - [`pybedtools`](https://daler.github.io/pybedtools/)
    - [`requests`](https://requests.readthedocs.io/en/master/)
    - [`scikit-learn`](https://scikit-learn.org/stable/install.html)
    - [`werkzeug`](https://werkzeug.palletsprojects.com/en/2.2.x/)

Additionally, to generate the example data, it requires:
* [`BigWigToBedGraph`](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
* [`MACS2`](https://github.com/macs3-project/MACS)

## Configuration

```
conda create -n ontarget -c bioconda -c conda-forge \
    biopython=1.79 \
    click=8.1.3 click-option-group=0.5.3 \
    flask=2.2.2 flask-cors=3.0.10 flask-restful=0.3.9 \
    fuzzywuzzy=0.18.0 \
    genomepy=0.14.0 \
    interval-binning=1.0.0 \
    nodejs=17.8.0 \
    numpy=1.23.5 \
    pandas=1.5.2 \
    pybedtools=0.9.0 \
    pymysql=1.0.2 \
    python=3.9.15 \
    python-levenshtein=0.20.9 \
    requests=2.28.1 \
    scikit-learn=1.2.0 \
    sqlalchemy=1.4.45 sqlalchemy-utils=0.38.3 \
    uwsgi=2.0.20 \
    werkzeug=2.2.2

pip install Flask-Limiter SQLAlchemy-FullText-Search
```

```python
import json
from GUD import GUDUtils

GUDUtils.user = "gud_r"
GUDUtils.pwd = ""
GUDUtils.host = "gud.cmmt.ubc.ca"
GUDUtils.port = 5506
GUDUtils.db = "mm10"
session = GUDUtils.get_session()

from ontarget.gene2interval import (get_intervals_limit_by_gene,
                                    get_intervals_limit_by_distance)

intervals = get_intervals_limit_by_gene(session, "Pitx3")
print(json.dumps(intervals, indent=4))
```
```
[
    {
        "chrom": "19",
        "start": 46135684,
        "end": 46152557,
        "type": "Gene interval",
        "id": "Pitx3 (from=Elovl3 to=Gbf1)",
        "score": 0,
        "strand": 0,
        "qualifiers": null
    }
]
```
```python
intervals = get_intervals_limit_by_distance(session, "Pitx3", 50000)
print(json.dumps(intervals, indent=4))
```
```
[
    {
        "chrom": "19",
        "start": 46085684,
        "end": 46198325,
        "type": "Gene interval",
        "id": "Pitx3 +/- 50 kb",
        "score": 0,
        "strand": 0,
        "qualifiers": null
    }
]
```
```python
from ontarget.interval2regions import get_regions

evidence = [
    ['./data/examples/Pitx3/dna_accessibility.bed.gz', 1.0],
    ['./data/examples/Pitx3/histone_modification-H3K27ac.bed.gz', 1.0],
    ['./data/examples/Pitx3/histone_modification-H3K36me3.bed.gz', 1.0],
    ['./data/examples/Pitx3/histone_modification-H3K4me1.bed.gz', 1.0],
    ['./data/examples/Pitx3/histone_modification-H3K4me3.bed.gz', 1.0],
    ['./data/examples/Pitx3/histone_modification-H3K9ac.bed.gz', 1.0]
]

regions = get_regions(session, chrom="19", start=46085684, end=46198325,
                      genome="mm10", evidence=evidence, liftover="hg19",
                      use_conservation=True, mask_exons=True,
                      mask_repeats=True)
print(json.dumps(regions, indent=4))
```
```
[
    ...
    {
        "chrom": "10",
        "start": 104000618,
        "end": 104001855,
        "type": "Promoter",
        "id": "RegulatoryRegion8",
        "score": 0.6321905313786568,
        "strand": 0,
        "qualifiers": {
            "source": "OnTarget",
            "genome": "mm10",
            "TSS": [
                [
                    "Pitx3",
                    -1
                ]
            ],
            "original coordinate": "mm10:chr19:46147762-46148962",
            "liftOver genome": "hg19",
            "sequence": "AGGCGCACCCCGGCTCCCAGCCCCAGGAGTCTGAGCCTAGGCCGAGGGGCATCGGGCGCCGTCAACTGCCCCCACACTGGGGGCCGCCCCGCCAGCTCCCCCGCTGGCCAGGGCCTCTGTAACCCCTTTCGCTCCCTGGCCGCCTGCCCCGTGCCCTGGGATTCCAGGGCGCTTTCTGGTAATGGTTTCCAGACTCCCTACTCCCTGCCTCGCGAAGCCCCTATCCTGGTTTCAAGTCTCCGGCTTCGTTCCGGATCCTCTGAGTCTCCCCCTTCCCTATCCCCTGGTTCTTACCCTGAGTTTCCCTTTTCCCAAGTCCTGGCCCCAGGTCTGCAGTCCGCCCTCCCATTCTCGACCTGTTCCCAAGGCACCGGAGTTCAGCTGCCCCAGGTCTCGTTTTAGGGATTCCAAGGTCCAGCAATAGCTCCTCGGCCCCATGAGCCCCTGTCCTTAGAATGGTCAAACACTCACAGTGCGTCCTGCAACAGGCAGACTCCCAGTAGCGGCGGCTGCGGCGGCGATCTAGAGGGCAGGCAGGGGCCAGGGGCCGGGCACCCGGCCGGAGTGGGGGCCGCCCCCCTGCTCCCGGGCCGCCTCTCCGCTCGGGCGCTCCTGGACTCTCGGAGGGAGTGAGCCTCACCGCGTACTGCCACCCCCAGCCGGCGCCCATTCACTTTATGGCAGACCAGGGCGCCCCCAGCCCGCCGCGGCGAGCCGCGCGCGTCAGGCCCCGCCCCTTTCCAGCTGCCCTGCTGGGGCTCCGCCCTTTCCAGCTGTGGATCTCCAGGCCCCGCCTTTGAGGGAGGGGTCTGGCCGGCGAGACGCCAAGAACCCCGCCCTCTGGCCAATCAGAAGCGCTCTTCAGCAACTCGGCCGCGCCCCTCCCCACGTGGCAGAGACCGCGCTCCGGCTAGGACGCTTAGGCAGAGCCAAGTGGGCGAGAGTAGAGTGGTCCCGGTAGCCACGGGTAAAAGGGATCGGCTGGCAGCGAAGTGGGTTGGGCCGCTACCCCGAGGAAAGTGGATCGGCGCCAGGTGGGAAGGCGCACGGGGCCGCCGGGTCTCAGCGCTCAGACTTCTCTGCCACTCAAGTTTCGCCGATTTTCTCCACTTGTCGCTGTCCGTCTCTCCTCCCTTTTGTCTCTGTTCCAATCCTCCCACGTCTCGCCTTTCTTTCCTCTTCTTTCTCCTCAGCCGCTTCGCTCCTGCTGTTTACCATTCTCCTTCCTCTGCCACTT",
            "enzymes": [
                "ABASI",
                "SAPI",
                "MAUBI",
                "CSP6I",
                ...
            ],
            "tfs": [
                "LBX1",
                "ZNF189",
                "FOXL1",
                "STAT6",
                ...
            ]
        }
    },
    ...
]
```
```python
from ontarget.regions2minips import get_minipromoters

minips = get_minipromoters(regions, enzymes=set(["AscI", "FseI"]),
                           tfs=set(["NR4A2", "PITX3"]))
print(json.dumps(minips, indent=4))
```
```
[
[
    {
        "chrom": "10",
        "starts": [
            103985959,
            103990858,
            103995880,
            104000233,
            104005794,
            104028873,
            104046030,
            104046360
        ],
        "ends": [
            103986305,
            103991344,
            103995998,
            104000416,
            104006152,
            104029024,
            104046180,
            104046597
        ],
        "type": "MiniPromoter",
        "ids": [
            "RegulatoryRegion1",
            "RegulatoryRegion4",
            "RegulatoryRegion6",
            "RegulatoryRegion7",
            "RegulatoryRegion11",
            "RegulatoryRegion15",
            "RegulatoryRegion16",
            "RegulatoryRegion17"
        ],
        "score": 0.5920698894940838,
        "strand": 1,
        "qualifiers": {
            "enzymes": [
                "ACCBSI",
                "SMAI",
                "ECORI",
                "HINDIII",
                ...
            ],
            "sequence": "CTTTAAATCCCCCCTCTTTTTCCTGTTGTAAAGGTCAGATTTTTTTTGACAACTGTGTGCAGGAACAAGATTTGCTAAATCCCTTCATTTTTATTCATACAAATATCGGAATTCAGTGAGCATCAGTTATTCCAGCATTTAGAGGTAAGGAAAACTCAGGGAGAATCTTGTGAATCTCAGCTTTGTGCTAGACATATATCTCTGGGGATAAATTTGCTGCCTCAACTTTAACCCCTTCATTCATTAAACCACTCCTCGCATTCCCTGTTTTCAAAGCTTATCAGGATGACTTTGATTTAGTTATTTGGAACCAAAAAGGCAGTTTGTGGAGCCTGGTCTGGAGCTGTGTAGCAGTAATATATCAGAGATCAGACTTAGGTCCCAGTGATGTCTCTTCTGTAACATCAGCATTCTGTTACCACTGTTACCACTGTAGGTCGAGTTGTGAAGGGCTGCCACATGTGCACAGATCACAACTGACCTTTCCTTGTTCCACGCTTATTGAGCAATTACTTGCCCGTTTGTGGAACTCACATTATCCCAGGCTTGACCTAGTATTTAGCTTCTCCTTTCCAGCGCTTTCTGCATCCTGAAGAGCCCAGTACTCAAAAGCTGAAAGCAAACCTAAGGAGCATTGTGTCTACTTTTATGGGTTTTCAGACTTTAGATGTGAAAGTCTCTCTAGATTCTGTACAGCGTCTACTTCTTTGTGATTGCCTCACTAATACAAAGGAACTGCTTCATTAATAATGATGAACTAGTGCACCCTAGAATCAGGTTGTCAGCAGAGGGACGTGTGGAAGGCTTGAGGTTATAGCAGCGCAGGGAACTGGGCTATGAGAAGGCTTCCCTCTTACAGGACAGCTCAGCTGCCAGTCTGGGTCTGTATGTAAACACTGGGATTTGAGGTGACCCACCGCAGCCGACCCCATCGCACACGCAGCAGCTGGGCAGGCGCTGCTGGGTAGACGGACATGCCACCCAGCCAGGGACATTCTGTCCGCCACTGCCAGACACAACTCCCGCTCCGGAGAGACCCACTTCGCTACCCGGACATTCTGTCCGCTCCGAGCCCCATGGGGCTGGTCAGGCAAAGGCTGGGGGCGCCAAGCTGTTGTCTGAGCTGGAGACAGGCTGGCTGAAGATTAAGCTCCCCCAAACTGCTGGAGGTCAGGAAGCTAGGGAGTTGGGACTGAGCTGCGGGCACGGGAGAAAGGCGGTCAGGGCCCGGGGCCGGGTCCCAGCGGCTGAAGGGCGGGCCACCCCGACGGGGCTTCCAGCCGGGGCCCTGCGGTCAATAAACGAGATGAGGTGGCTAGAGACGGGGTGGGAGGGAAGGAGCGCACGGAGTCCGGGGTCGGAGAACAGAAGGCGCGGGCGGCTGGACCGAGTCATCGGCGCCAGGGTCTGGCGGGACAGTAGGATGGGGTTGAGGAGGGGTGGTGAGGGGAAGGAGAGACGGTGTCAGGGCCGAAGAAGCGCGCGGTCCGAGTAGTAGGGAGCAGTGGGAGCAGAGGCTGGAGGTTGGCGGGCAAGAGACCGCGTGGGGCTGCGGCCGGGAGCCAGGGTCCGGGGTCCGGGGTCCGAGGGAGGGGGCAGGTGGGGTGGAACCGCTGGCCTCCGGGTCGCAGGCTGAGCGCGGAGGGCCCGCGCGGGTGCGAGTCGCGGGTCTGGAGAGCATACGCAGGCGGGCCCTTGAGTATGACCAATCAGAATGCGGACTGCGTCCCAGGGGCGGAGCAGAGGCGTATCTTGGTCGAGATTGGATAGCGGCGGGGCGCAGGAAAGAGGTCGCGCCAGCCCGGGCAGGCAGCTTTGCAAGTCCGCGTTATATATCGCAGTGGCTGCGCCCGGGATAGCTGGCTGCGCCGCCGCGCACATGCCTAGGTTCGACGCCCTCCTCCCTTTGCCCAGGAGTTCCTTCTGTCCCGGCTCTGTTCCGTCTCGCCCCGAGGTTCACGCCATCCTCGGAGCCCCAGCCTTTCACCCAGCGCCTCCAAGCTTTGGACCTTGACTTCTGCAAAACTAG",
            "size": 2029,
            "tfs": [
                "AHR",
                "ALX3",
                "AR",
                "ARGFX",
                ...
            ]
        }
    },
    ...
]
```
