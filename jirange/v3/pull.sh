#!/bin/bash

set -ex

for fn in 20230325_k21_008kbit_out.txt \
	  20230325_k21_032kbit_out.txt \
	  20230325_k21_128kbit_out.txt \
	  20230325_k31_008kbit_out.txt \
	  20230325_k31_032kbit_out.txt \
	  20230325_k31_128kbit_out.txt
do
    scp blangme2@login.rockfish.jhu.edu:${fn} .
done
