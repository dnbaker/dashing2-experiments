#!/bin/bash

#
# Download all the pairs of genomes needed for the JI-range experiments
#
# Put them all, still compressed, in subdirectories like pairs/0001/1.fna.gz /
# pairs/0001/2.fna.gz
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
        dir="pairs/${npad}"
        mkdir -p "${dir}"
        pushd "${dir}"
        fn="${i}.fna.gz"
        URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/${gcf}/${tri1}/${tri2}/${tri3}/${i_wg}/${fn}"
        if [[ -f "${fn}" ]] ; then
            mv "${fn}" "${field}.fna.gz"
        elif [[ ! -f "${field}.fna.gz" ]] ; then
            wget "${URL}" -O "${field}.fna.gz"
        fi
        popd
        n=$(expr $n + 1)
    done
done

