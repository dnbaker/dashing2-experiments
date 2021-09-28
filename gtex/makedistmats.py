from glob import glob
from multiprocessing import cpu_count as CC
import numpy as np
import sketch



ss = glob("GTexSplicing.stackedsketches.*indices*i64")

pc = sketch.util.pcount_eq

mms = list(map(lambda x: np.memmap(x, dtype=np.uint64).reshape(19214, -1), ss))


NT = CC()

def make_distmats():
    ret = []
    for mm, sf in zip(mms, ss):
        distmat = sketch.util.pcount_eq(mm, NT)
        df = sf + ".distmat.u32"
        distmat.tofile(df)
        ret.append(df)
    return ret

dms = make_distmats()
