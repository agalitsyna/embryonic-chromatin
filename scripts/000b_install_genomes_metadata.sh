cd ./data/genomes/
genomepy plugin enable bwa

### Manually curated files provided with the repo:
# danRer11.armsizes.txt
# danRer11.centromeres.manual.txt
# danRer11.reduced.chromsizes

### Danio rerio genome:
# genome:
wget http://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz -P danRer11
# filter out unconventional chromosomes:
gzip -d danRer11/danRer11.fa.gz
samtools faidx danRer11/danRer11.fa
bedtools getfasta -nameOnly -fi danRer11/danRer11.fa -bed <(awk '{print $1, "0", $2, $1}' danRer11/danRer11.reduced.chromsizes | sed 's/ /\t/g') > danRer11/danRer11.reduced.fa
# build bwa index:
bwa index danRer11/danRer11.reduced.fa

### Mouse:
genomepy install mm10 -p UCSC -g ./

### Human:
genomepy install hg38 -p UCSC -g ./

### Medaka:
genomepy install oryLat2 -p UCSC -g ./

### Xenopus tropicalis:
genomepy install xenTro10 -p UCSC -g ./

### Drosophila melanogaster:
genomepy install dm6 -p UCSC -g ./