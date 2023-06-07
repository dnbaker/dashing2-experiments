#!/bin/bash

#
# Don't forget to 'singularity run docker://benlangmead/dashing2-experiments'
#

set -ex

# All sketches will be 8K, i.e. 1000 8-byte estimators

NJOBS=16
K=21
SIZE=1024

#
# Sourmash sketching
#
if [[ ! -f sourmash_sketch.time ]] ; then
  python -c "for i in range(1010): print('%04d' % (i+1))" | \
    /usr/bin/time -v parallel --jobs ${NJOBS} \
      sourmash sketch dna -p k=${K},noabund,num=${SIZE} {}_1.fna.gz {}_2.fna.gz -o {}_k${K}.sig 2>&1 | \
      tee sourmash_sketch.time
fi

#
# Mash sketching
#
if [[ ! -f mash_sketch.time ]] ; then
  /usr/bin/time -v \
    mash sketch -s ${SIZE} -p ${NJOBS} -k ${K} -o all.msh *.fna.gz 2>&1 | \
    tee mash_sketch.time
fi

#
# Dashing 2 OP sketching
#
if [[ ! -f d2_sketch.time ]] ; then
  ls *.fna.gz > fastas.txt
  /usr/bin/time -v \
    dashing2 sketch --cache -k${K} --oneperm -S${SIZE} -p${NJOBS} -F fastas.txt 2>&1 | \
    tee d2_sketch.time
fi

#
# Dashing 2 full sketching
#
if [[ ! -f d2f_sketch.time ]] ; then
  ls *.fna.gz > fastas.txt
  /usr/bin/time -v \
    dashing2 sketch --cache -k${K} --full-setsketch -S${SIZE} -p${NJOBS} -F fastas.txt 2>&1 | \
    tee d2f_sketch.time
fi

#
# Dashing 2 weighted sketching
#
if [[ ! -f d2w_sketch.time ]] ; then
  ls *.fna.gz > fastas.txt
  /usr/bin/time -v \
    dashing2 sketch --cache --prob --countsketch-size 5000000 -S${SIZE} -p${NJOBS} -F fastas.txt 2>&1 | \
    tee d2w_sketch.time
fi

#
# Sourmash all-pairs Jaccard
#
if [[ ! -f sourmash_cmp.time ]] ; then
  /usr/bin/time -v \
    sourmash compare *.sig -o sm.np -p ${NJOBS} 2>&1 | \
    tee sourmash_cmp.time
fi

#
# Mash all-pairs Jaccard
#
if [[ ! -f mash_cmp.time ]] ; then
  /usr/bin/time -v \
    mash triangle -p ${NJOBS} all.msh 2>&1 | \
    tee mash_cmp.time
fi

#
# Dashing 2 OP all-pairs Jaccard
#
if [[ ! -f d2_cmp.time ]] ; then
  #nsketches=$(ls *.opss *.ss *.pmh 2>/dev/null | wc -l)
  ls *.opss > opss.txt
  /usr/bin/time -v \
    dashing2 cmp -p${NJOBS} --presketched -F opss.txt --cmpout d2_cmp.fl32 2>&1 |
    tee d2_cmp.time
  #nsketches_new=$(ls *.opss *.ss *.pmh 2>/dev/null | wc -l)
  #if [[ ${nsketches} != ${nsketches_new} ]] ; then
  #  echo "Oops, created more sketches with cmp command"
  #fi
fi

#
# Dashing 2 full all-pairs Jaccard
#
if [[ ! -f d2f_cmp.time ]] ; then
  ls *.ss > ss.txt
  /usr/bin/time -v \
    dashing2 cmp -p${NJOBS} --presketched -F ss.txt --cmpout d2f_cmp.fl32 2>&1 |
    tee d2f_cmp.time
fi

#
# Dashing 2 weighted all-pairs Jaccard
#
if [[ ! -f d2w_cmp.time ]] ; then
  ls *.pmh > pmh.txt
  /usr/bin/time -v \
    dashing2 cmp -p${NJOBS} --presketched -F pmh.txt --cmpout d2w_cmp.fl32 2>&1 |
    tee d2w_cmp.time
fi
