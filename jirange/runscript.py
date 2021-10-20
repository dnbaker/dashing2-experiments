import sys
import subprocess
from argparse import ArgumentParser as AP

ap = AP()
ap.add_argument("input_files", nargs="+", help="Paths to concatenate into a table for analysis and plotting")

args = ap.parse_args()

subprocess.check_call(["Rscript", "00wide-large.R"] + args.input_files)
for script in '01pivoted-large.R,02augmented-large.R,03filtered-large.R,04sse-large.R,06augmented_subset-large.R,07filtered_subset-large.R,08sse_subset-large.R,10anisumm-large.R,20sedf-large.R,30boxplot-large.R,40pairscatter-large.R'.split(','):
    subprocess.check_call(["Rscript", script])
