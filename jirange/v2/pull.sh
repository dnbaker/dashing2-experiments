#!/bin/bash

set -ex

for fn in 20220928_k21_008kbit_out.txt \
	      20220928_k21_032kbit_out.txt \
	      20220928_k21_128kbit_out.txt \
	      20220928_k31_008kbit_out.txt \
	      20220928_k31_032kbit_out.txt \
	      20220928_k31_128kbit_out.txt
do
    scp blangme2@login.rockfish.jhu.edu:/home/blangme2/scr16_blangme2/langmead/dashing/${fn} .
done
