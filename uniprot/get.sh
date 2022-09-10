#!/bin/bash

set -ex

if [[ ! -f uniprot_sprot.fasta ]] ; then
  wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
  gzip -dc uniprot_sprot.fasta.gz > uniprot_sprot.fasta
fi

# 568002 for 2022_03
echo "Uniprot:"
grep -c '^>' uniprot_sprot.fasta


if [[ ! -f uniref50.fasta ]] ; then
  wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
  gzip -dc uniref50.fasta.gz > uniref50.fasta
fi

# 53625855 for 2022_03
echo "Uniref 50:"
grep -c '^>' uniref50.fasta
