import sys
from collections import defaultdict
import numpy as np


def lines2np(x):
    return np.array(list(map(float, map(lambda x: x.split('\t')[3], x))))

if __name__ == "__main__":
    lines = list(filter(lambda x: not x.isspace(), open(sys.argv[1])))
    hlls, ssd, ssh, ssb, ssn = lines[::5], lines[1::5], lines[2::5], lines[3::5], lines[4::5]
    hest = []
    for x in hlls:
        try:
            hest.append(float(x.split('\t')[3]))
        except:
            print(x)
            raise
    hest, csde, cshe, csbe, csne = map(lines2np, (hlls, ssd, ssh, ssb, ssn))
    tv = np.array([float(x.split('\t')[4]) for x in ssb])
    fullarr = np.vstack([hest, csde, cshe, csbe, csne, tv])
    diffs = fullarr[:,:].T - tv[:,np.newaxis]
    absdiffs = np.abs(diffs)
    reldiffs = diffs / tv[:,np.newaxis]
    absreldiffs = absdiffs / tv[:,np.newaxis]
    mae = np.mean(absdiffs, axis=0)
    mede = np.median(absdiffs, axis=0)
    mede2 = np.median(np.square(absdiffs), axis=0)
    mb = np.mean(diffs, axis=0)
    ard = np.mean(absreldiffs, axis=0)
    ard2 = np.mean(np.square(absreldiffs), axis=0)
    medrd = np.median(absreldiffs, axis=0)
    medrd2 = np.median(np.square(absreldiffs), axis=0)
    print("#Measure\tHLL\tCSD\tCS16\tCS8\tCS4")
    print("\t".join(["MeanAbsEr"] + list(map(str, mae[:5]))))
    print("\t".join(["MeanBias"] + list(map(str, mb[:5]))))
    print("\t".join(["MeanRelDiff"] + list(map(str, ard[:5]))))
    print("\t".join(["MeanRelDiffSquared"] + list(map(str, ard2[:5]))))
    print("\t".join(["MedianRelDiff"] + list(map(str, mede[:5] / np.max(mede)))))
    print("\t".join(["MedianRelDiffSquared"] + list(map(str, mede2[:5] /np.max(mede2)))))
    print("\t".join(["MeanAbsErr_vsHLL"] + list(map(str, mae[:5] / np.max(mae)))))
    print("\t".join(["MeanBias_vsHLL"] + list(map(str, mb[:5] / np.max(mb)))))
    print("\t".join(["MeanRelDiff_vsHLL"] + list(map(str, ard[:5] / np.max(ard)))))
    print("\t".join(["MeanRelDiffSquared_vsHLL"] + list(map(str, ard2[:5] /np.max(ard2)))))
    print("\t".join(["MedianRelDiff_vsHLL"] + list(map(str, medrd[:5] / np.max(medrd)))))
    print("\t".join(["MedianRelDiffSquared_vsHLL"] + list(map(str, medrd2[:5] /np.max(medrd2)))))
