#!/bin/bash

set -ex

if [[ ! -f sourmash.time ]] ; then
  python -c "for i in range(1010): print('%04d' % (i+1))" | \
    /usr/bin/time -v parallel --jobs 10 \
      sourmash sketch dna -p k=${k},noabund,num=1000 {}_1.fna.gz {}_2.fna.gz -o {}_k${k}.sig 2>&1 | \
      tee sourmash.time
fi
