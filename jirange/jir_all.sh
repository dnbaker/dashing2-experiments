#!/bin/bash

set -ex

JIR_SCR="dashing2-experiments/jirange/jir.py"

test -f "${JIR_SCR}"

export PATH=$PATH:$HOME/software

DATE="20230325"

# 21-mers, 8 kbits
name="${DATE}_k21_008kbit"
python3 "${JIR_SCR}" bva_buckets.tsv bva_fastas.txt -k21 -s13 -S14 --name ${name} -o ${name}_out.txt--cpu 20
# 21-mers, 32 kbits
name="${DATE}_k21_032kbit"
python3 "${JIR_SCR}" bva_buckets.tsv bva_fastas.txt -k21 -s15 -S16 --name ${name} -o ${name}_out.txt --cpu 20
# 21-mers, 128 kbits
name="${DATE}_k21_128kbit"
python3 "${JIR_SCR}" bva_buckets.tsv bva_fastas.txt -k21 -s17 -S18 --name ${name} -o ${name}_out.txt --cpu 20

# 31-mers, 8 kbits
name="${DATE}_k31_008kbit"
python3 "${JIR_SCR}" bva_buckets.tsv bva_fastas.txt -k31 -s13 -S14 --name ${name} -o ${name}_out.txt --cpu 20
# 31-mers, 32 kbits
name="${DATE}_k31_032kbit"
python3 "${JIR_SCR}" bva_buckets.tsv bva_fastas.txt -k31 -s15 -S16 --name ${name} -o ${name}_out.txt --cpu 20
# 31-mers, 128 kbits
name="${DATE}_k31_128kbit"
python3 "${JIR_SCR}" bva_buckets.tsv bva_fastas.txt -k31 -s17 -S18 --name ${name} -o ${name}_out.txt --cpu 20
