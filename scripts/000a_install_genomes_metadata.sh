cd ./data/genomes/
genomepy plugin disable bwa || true

### Manually curated files provided with the repo:
# danRer11.armsizes.txt
# danRer11.centromeres.manual.txt
# danRer11.reduced.chromsizes

### Danio rerio genome:
# genes:
wget http://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/genes/danRer11.ensGene.gtf.gz -P danRer11


