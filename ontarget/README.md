```
./gene2interval.py Adora2a mm10 --lim-by-dist 100 -o ../data/examples/Adora2a/
./interval2regions.py --conservation --mask-exons --mask-repeats \
    --evidence ../data/examples/Adora2a/dna_accessibility.bed.gz 1. \
    --evidence ../data/examples/Adora2a/histone_modification-H3K27ac.bed.gz 1. \
    --evidence ../data/examples/Adora2a/histone_modification-H3K36me3.bed.gz 1. \
    --evidence ../data/examples/Adora2a/histone_modification-H3K4me1.bed.gz 1. \
    --evidence ../data/examples/Adora2a/histone_modification-H3K4me3.bed.gz 1. \
    --evidence ../data/examples/Adora2a/tf_binding-POL2.bed.gz 1. \
    --liftover hg19 -o ../data/examples/Adora2a/ 10 75216876 75434792 mm10
./regions2minips.py --enzyme AscI --enzyme FseI --tf RARB --tf SP9 \
    -o ../data/examples/Adora2a/ ../data/examples/Adora2a/regions.json

./gene2interval.py Pitx3 mm10 --lim-by-dist 50 -o ../data/examples/Pitx3/
./interval2regions.py --conservation --mask-exons --mask-repeats \
    --evidence ../data/examples/Pitx3/dna_accessibility.bed.gz 1. \
    --evidence ../data/examples/Pitx3/histone_modification-H3K27ac.bed.gz 1. \
    --evidence ../data/examples/Pitx3/histone_modification-H3K36me3.bed.gz 1. \
    --evidence ../data/examples/Pitx3/histone_modification-H3K4me1.bed.gz 1. \
    --evidence ../data/examples/Pitx3/histone_modification-H3K4me3.bed.gz 1. \
    --evidence ../data/examples/Pitx3/histone_modification-H3K9ac.bed.gz 1.\
    --liftover hg19 -o ../data/examples/Pitx3/ 19 46085684 46198325 mm10
./regions2minips.py --enzyme AscI --enzyme FseI --tf NR4A2 --tf PITX3 \
    -o ../data/examples/Pitx3/ ../data/examples/Pitx3/regions.json
```