#!/bin/bash

set -ex

DASH2="${HOME}/scr16_blangme2/langmead/dashing2/dashing2-64"

test -f uniprot_sprot.fasta

test -x "${DASH2}"

# Ran on devlangmead1
/usr/bin/time -v \
  "${DASH2}" \
  sketch --parse-by-seq --topk 256 -S256 -p80 --cmpout out.txt -k 5 --protein6 uniprot_sprot.fasta
