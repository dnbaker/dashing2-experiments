import sys
from collections import defaultdict
import numpy as np

if __name__ == "__main__":
    lines = list(open(sys.argv[1]))
    hlls, ssd, ssh, ssb = lines[::4], lines[1::4], lines[2::4], lines[3::4]
    hest = np.array([float(x.split('\t')[3]) for x in hlls])
    csde = np.array([float(x.split('\t')[3]) for x in ssd])
    cshe = np.array([float(x.split('\t')[3]) for x in ssh])
    csbe = np.array([float(x.split('\t')[3]) for x in ssb])
    tv = np.array([float(x.split('\t')[4]) for x in ssb])
    fullarr = np.vstack([hest, csde, cshe, csbe, tv])
    diffs = fullarr[:,:].T - tv[:,np.newaxis]
    absdiffs = np.abs(diffs)
    reldiffs = diffs / tv[:,np.newaxis]
    absreldiffs = absdiffs / tv[:,np.newaxis]
    mae = np.mean(absdiffs, axis=0)
    mb = np.mean(diffs, axis=0)
    ard = np.mean(absreldiffs, axis=0)
    ard2 = np.mean(np.square(absreldiffs), axis=0)
    print("#Measure\tHLL\tCSD\tCS16\tCS8")
    print("\t".join(["MeanAbsEr"] + list(map(str, mae[:4]))))
    print("\t".join(["MeanBias"] + list(map(str, mb[:4]))))
    print("\t".join(["MeanRelDiff"] + list(map(str, ard[:4]))))
    print("\t".join(["MeanRelDiffSquared"] + list(map(str, ard2[:4]))))
    print("\t".join(["MeanAbsErr_vsHLL"] + list(map(str, mae[:4] / np.max(mae)))))
    print("\t".join(["MeanBias_vsHLL"] + list(map(str, mb[:4] / np.max(mb)))))
    print("\t".join(["MeanRelDiff_vsHLL"] + list(map(str, ard[:4] / np.max(ard)))))
    print("\t".join(["MeanRelDiffSquared_vsHLL"] + list(map(str, ard2[:4] /np.max(ard2)))))
