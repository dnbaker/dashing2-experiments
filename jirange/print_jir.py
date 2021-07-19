from sys import argv, stdout, stderr
import numpy as np
import argparse
# header = "#ANI\tWJI\tJI\tMash\tDash1\tBD8\tBD4\tBD2\tBD1\tBDN\tSS8\tSS2\tSS1\tSSN\tFSS8\tFSS2\tFSS1\tFSSN\tMH8\tMH4\tMH2\tMH1\tMHN"
# Newest version, including full and one-perm
header = "#ANI\tWJI\tJI\tMash\tDash1\tBD8\tBD4\tBD2\tBD1\tBDN\tSS8\tSS2\tSS1\tSSN\tFSS8\tFSS2\tFSS1\tFSSN\tMH8\tMH4\tMH2\tMH1\tMHN\tFMH8\tFMH4\tFMH2\tFMH1\tFMHN"
ap = argparse.ArgumentParser()
ap.add_argument("matrix", help="Emitted numpy array from jirange/jir.py, which is (x by 28) in shape")
ap.add_argument("settings", help="Emitted results from jirange/jir.py, which has tuples for the experimental conditions")
ap = ap.parse_args()

mat = np.memmap(ap.matrix, np.float32).reshape(-1, 23)

print(header)
for mr, l in zip(mat, open(ap.settings)):
    lh, rh, k, ss = l.strip().split('\t')
    fstr = "\t".join(map(str, mr))
    print(f"{lh}\t{rh}\t{k}\t{ss}\t{fstr}", file=stdout)
