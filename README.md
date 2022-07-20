# k-mer based phylogenetic reconstruction

A simple proof-of-concept of reconstructing history from publicly
available reference genomes.

## Cloning the directory

```
git clone https://github.com/guilhermesena1/phylo.git
cd phylo
```

## Compiling the code
```
make
```

## Downloading genomes

```
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/bigZips/gorGor6.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/criGriChoV2/bigZips/criGriChoV2.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/felCat9/bigZips/felCat9.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/vicPac2/bigZips/vicPac2.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/balAcu1/bigZips/balAcu1.fa.gz
```

```
for i in $(cat genome_urls.txt); do echo $i; wget ${i}; done
gunzip *.fa.gz;
mkdir genomes
mv *.fa genomes
```

## Running the reconstruction in C++
```
./phylo genome_inputs.txt >kmer-counts.tsv
```

## Creating the hierarchical clustering in R
```
> x <- read.table('kmer-counts.tsv', header = T, row.names=NULL)
> plot(hclust(dist(t(x))), hang = -1, xlab = "species", ylab = "k-mer squared distance")
```
