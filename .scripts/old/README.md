pip install biopython
pip install mysqlclient
pip install pybedtools
pip install sqlalchemy


#Th design
/usr/local/python2.7.14/bin/python ./getRR.py -gs Th -c midbrain -o ./genes/Th/Th.bed -v
#MSN designs
/usr/local/python2.7.14/bin/python ./getRR.py -gs GENE_SYMBOL -f midbrain -o ./genes/GENE_SYMBOL/GENE_SYMBOL.bed -v
