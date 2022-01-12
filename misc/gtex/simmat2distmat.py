from glob import glob
import numpy as np


dms = glob("*distmat.u16")
assert len(dms) == 1

def sz2t(x):
    if x < 256:
        return np.uint8
    if x < 65536:
        return np.uint16
    if x < (1<<32):
        return np.uint32
    return np.uint64

sketchsizes = [int(x.split(".")[-4]) for x in dms]
dtypes = list(map(sz2t, sketchsizes))
sketchinvs = [1. / x for x in sketchsizes]

for path, ssz, dt, siv in zip(dms, sketchsizes, dtypes, sketchinvs):
    dm = np.memmap(path, dtype=dt)
    distmat = np.zeros(dm.shape, dtype=np.float32)
    nz = np.where(dm)
    zeros = np.where(np.logical_not(dm))
    distmat[zeros] = np.inf
    distmat[nz] = -np.log(dm[nz] * siv)
    distmat.tofile(path[:-3] + "f32")
