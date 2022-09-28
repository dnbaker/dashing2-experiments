import numpy as np
import numba
import argparse as AP
from sys import stderr


ap = AP.ArgumentParser()
ap.add_argument("distmat", help="Path to a binary output of pairwise Jaccard similarities, in float32")
ap.add_argument("--nbuckets", "-n", type=int, help="Number of buckets to divide the data into. Default: 100.", default=100)
ap.add_argument("--fillnum", "-F", type=int, help="Default: 10.", default=10)
ap = ap.parse_args()

nb = ap.nbuckets
mat = np.memmap(ap.distmat, np.float32)
nelem = int(np.ceil(np.sqrt(np.prod(mat.shape) * 2)))


@numba.cfunc("uint32(double)")
def id2bkt(x):
    if x == 0.: return 0
    return 1 + int(x * nb)


def maketi(n):
    lh, rh = np.triu_indices(n, 1)
    return list(zip(lh, rh))


id2bkt = np.vectorize(id2bkt)

indices = maketi(nelem)

buckets = [[] for i in range(nb + 2)]

seen = set()
mat_num_registers = len(mat)
nbfull = 0
fillnum = ap.fillnum
for ci, cv in enumerate(map(id2bkt, mat)):
    myb = buckets[cv]
    if len(myb) >= fillnum:
        continue
    x, y = indices[ci]
    if x in seen or y in seen:
        continue
    seen.add(x)
    seen.add(y)
    myb.append(indices[ci])  # maybe check it's not redundant?
    print("  appended: " + str(indices[ci]) + " to bucket " + str(cv), file=stderr)
    if len(myb) == fillnum:
        nbfull += 1
        print(f"{nbfull}/{len(buckets)}", file=stderr)
        if nbfull > len(buckets) - 10:
            remaining = []
            for i, bucket in enumerate(buckets):
                if len(bucket) < fillnum:
                    remaining.append(str(i))
                else:
                    assert len(bucket) == fillnum
            print("  remaining buckets: " + ' '.join(remaining), file=stderr)
        if nbfull == len(buckets):
            print(f"Finished all, breaking after {ci}/{mat_num_registers} {ci * 100./mat_num_registers}", file=stderr)
            break

rstrs = ["0."]
print("Bucket 0 [0.0]\t" + "\t".join(map(lambda x: "%d:%d" % (x[0], x[1]), buckets[0])))


def idx2v(idx):
    return idx * (1. / nb), (idx + 1) * (1. / nb)


for idx, bk in enumerate(buckets[1:-1]):
    iv = idx + 1
    s, t = idx2v(idx)
    print(f"Bucket {iv} [{s}->{t}]\t" + "\t".join(map(lambda x: "%d:%d" % (x[0], x[1]), bk)))
print("Bucket " + str(nb) + " [1.]\t" + "\t".join(map(lambda x: "%d:%d" % (x[0], x[1]), buckets[-1])))
