from sys import argv, stdout, stderr
import numpy as np
import argparse
header = "#LHGenome\tRHGenome\tk\tSketchsize\tMash\tDashing1\tTrueJI\tTrueWJI\tSS8Bytes\tSS2Bytes\tSS1Byte\tSSNibble\tMH8Bytes\tMH4Bytes\tMH2Bytes\tMH1Byte\tMHNibble"
ap = argparse.ArgumentParser()
ap.add_argument("matrix", help="Emitted numpy array from jirange/jir.py, which is (x by 13) in shape")
ap.add_argument("settings", help="Emitted results from jirange/jir.py, which has tuples for the experimental conditions")
ap = ap.parse_args()

mat = np.memmap(ap.matrix, np.float32).reshape(-1, 13)

print(header)
for mr, l in zip(mat, open(ap.settings)):
    lh, rh, k, ss = l.strip().split('\t')
    fstr = "\t".join(map(str, mr))
    print(f"{lh}\t{rh}\t{k}\t{ss}\t{fstr}", file=stdout)
