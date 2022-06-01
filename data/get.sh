#!/bin/bash

set -ex
which wget

while IFS= read -r line
do
  tri1=${line:4:3}
  tri2=${line:7:3}
  tri3=${line:10:3}
  url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${tri1}/${tri2}/${tri3}/${line}/${line}_genomic.fna.gz"
  wget ${url}
done < manifest.csv
