This folder contains scripts for selecting `GUD` features from which to identify *cis*-regulatory regions (CRRs). As use example, we focus on how to obtain CRRs in `endothelial cells`.

### 1) Select relevant samples

Table `samples` enables fulltext searches in `IN NATURAL LANGUAGE MODE`.

```
python -m GUD.scripts.name2sample --name "endothelial cells" -o ./examples/samples.txt
```

or

```
from GUD import GUDglobals
from GUD.ORM.sample import Sample

# Establish a GUD session
session = GUDglobals.establish_GUD_session()

# Get samples
samples = Sample.select_by_fulltext(
    session,
    "endothelial cells"
)
```

Output samples are sorted by relevance and, in some cases, include not relevant samples (*e.g.* `endothelial progenitor cell (derived from CD14-positive monocyte)`).

### 2) Select differentially expressed genes

The script `sample2gene.py` identifies the transcriptions start sites (TSSs) of genes that are sufficiently expressed in the input samples (by default 100 tags per million; option `--tpm`) and whose expression levels in the input samples exceeds a minimum percentile over the total.

```
python -m GUD.scripts.sample2gene --sample-file ./examples/samples.txt -o  ./examples/genes.txt
```

or

```
from GUD import GUDglobals
from GUD.ORM.expression import Expression
from GUD.ORM.gene import Gene
from GUD.ORM.sample import Sample
from GUD.ORM.tss import TSS

# Establish a GUD session
session = GUDglobals.establish_GUD_session()

# Get samples
samples = Sample.select_by_fulltext(
    session,
    [
        'endothelial cell (aorta)',
        'endothelial cell (vein)',
        'endothelial cell (artery)',
        'endothelial cell (umbilical vein)'
    ]
)
```

Output TSSs (by default limited to a maximum of 50; option `--tss`) are sorted by their differential expression in the input samples.

### 2) Select TSS from samples
```
from GUD import GUDglobals
from GUD.scripts.sample2gene import get_differentially_expressed_tss

# Establish a GUD session
session = GUDglobals.establish_GUD_session()

# Get differentially expressed TSSs
diff_exp_tss = get_differentially_expressed_tss(
    session,
    samples=["GM12878", "B cell (CD19-positive)"]
)
```    

### 3) Select a region for a gene:
```
from GUD import GUDglobals
from GUD.scripts.gene2region import get_gene_region

# Establish a GUD session
session = GUDglobals.establish_GUD_session()

# Get that gene's region
region = get_gene_region(
    session,
    gene="RELA",
    samples=["GM12878", "B cell (CD19-positive)"]
)
```