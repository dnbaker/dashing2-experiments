from sys import argv, stdout, stderr
import numpy as np
header = "#LHGenome\tRHGenome\tk\tSketchsize\tMash\tDashing1\tTrueJI\tTrueWJI\tSS8Bytes\tSS2Bytes\tSS1Byte\tSSNibble\tMH8Bytes\tMH4Bytes\tMH2Bytes\tMH1Byte\tMHNibble"

mat = np.memmap(argv[1], np.float32).reshape(-1, 13)

print(header)
for mr, l in zip(mat, open(argv[2])):
    lh, rh, k, ss = l.strip().split('\t')
    fstr = "\t".join(map(str, mr))
    print(f"{lh}\t{rh}\t{k}\t{ss}\t{fstr}", file=stdout)
