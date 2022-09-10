#!/bin/bash

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gzip -dc uniprot_sprot.fasta.gz > uniprot_sprot.fasta

# 568002 for 2022_03
grep -c '^>' uniprot_sprot.fasta