#!/bin/bash

#
# Download all the pairs of genomes needed for the JI-range experiments
#

set -ex

if [[ ! -f pairs1010.csv ]] ; then
    wget https://www.cs.jhu.edu/~langmea/resources/d2/pairs1010.csv
fi

for field in 2 1 ; do
    n=1
    for i in $(cut -d ',' -f ${field} pairs1010.csv) ; do
        i_wg=${i/%_genomic/}
        gcf=${i:0:3}
        tri1=${i:4:3}
        tri2=${i:7:3}
        tri3=${i:10:3}
        npad=$(printf %04d $n)
        fn="${i}.fna.gz"
        URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/${gcf}/${tri1}/${tri2}/${tri3}/${i_wg}/${fn}"
        wget "${URL}" -O "${npad}_${field}.fna.gz"
        n=$(expr $n + 1)
    done
done

