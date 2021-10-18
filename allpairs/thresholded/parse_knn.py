import os
import sys
import numpy as np


def parse_knn(path):
    ifp = open(path)
    names = [x.split()[0] for x in ifp][1:]
    name2id = {name: idx for idx, name in enumerate(names)}
    ifp.seek(os.SEEK_SET)
    next(ifp)
    return [np.fromiter(map(name2id.__getitem__, (x.split(":")[0] for x in line.split()[1:])), np.uint32) for line in ifp]


if __name__ == "__main__":
    out = parse_knn(sys.argv[1])
    print(out)
