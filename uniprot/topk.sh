#!/bin/bash

# ./topk.sh 2>&1 | tee topk.log

set -ex

DASH2="${HOME}/scr16_blangme2/langmead/dashing2/dashing2-64"

test -f uniprot_sprot.fasta

test -x "${DASH2}"

# Ran on devlangmead1

#
# UNIPROT
#

# This is an attempt to recreate the experiment that generates the first KNN
# graph mentioned in "All-pairs comparisons using LSH" section of Dashing 2
# paper.

/usr/bin/time -v \
  "${DASH2}" \
  sketch --parse-by-seq --topk 256 -S256 -p80 --cmpout uniprot_knn.tsv -k 5 --protein6 uniprot_sprot.fasta

# Now need experiment that generates the graph using full all-pairs experiment.


#
# UNIREF 50
#

/usr/bin/time -v \
  "${DASH2}" \
  sketch --parse-by-seq --topk 256 -S256 -p80 --cmpout uniref50_knn.tsv -k 5 --protein6 uniref50.fasta

# Now need experiment that generates the graph using full all-pairs experiment.
# (This one has to get interrupted after a while)

