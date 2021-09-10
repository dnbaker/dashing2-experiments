import sys
import os
import numpy as np
import subprocess
from time import time


def bch_sketch_mash(pathf, k, threads, destp, size):
    startt = time()
    subprocess.check_call(f"mash sketch -k {k} -o {destp} -l {pathf} -p {threads}", shell=True)
    stopt = time()
    return destp, (stopt - startt)


def bch_sketch_dashing(pathf, k, threads, size):
    startt = time()
    subprocess.check_call(f"dashing sketch -S {int(np.ceil(np.log2(size)))} -k {k} -F {pathf} -p {threads}", shell=True)
    stopt = time()
    return pathf, (stopt - startt)


def bch_sketch_dashing2(pathf, k, threads, size, oneperm=False):
    startt = time()
    pstr = " --oph " if oneperm else ""
    subprocess.check_call(f"dashing2 sketch -k {k} -S {size} -F {pathf} -p {threads} {pstr}", shell=True)
    return pathf, (time() - startt)


def bch_dist_dashing2(pathf, k, threads, size, oneperm=False, distdest=None, binary=False, regsize=-1):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    if regsize <= 0. or regsize not in [0.5, 1, 2, 4, 8]:
        raise ValueError("regsize must be > 0 and in [.5, 1, 2, 4, 8]")
    startt = time()
    pstr = " --oph " if oneperm else ""
    pstr += " --binary-output " if binary else ""
    subprocess.check_call(f"dashing2 sketch -k {k} --cmpout {distdest} -S {size} -F {pathf} -p {threads} {pstr}", shell=True)
    return pathf, (time() - startt)


def bch_dist_dashing(pathf, k, threads, size, distdest=None, binary=False):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    startt = time()
    pstr = " --emit-binary" if binary else ""
    subprocess.check_call(f"dashing dist -k {k} -O{distdest} -o {distdest + '.sizes'} -S {int(np.ceil(np.log2(size)))} -F {pathf} -p {threads} {pstr}", shell=True)
    return pathf, (time() - startt)


def bch_dist_bindash(pathf, threads, size, distdest=None):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    startt = time()
    subprocess.check_call(f"bindash dist --outfname={distdest} --nthreads={threads} {pathf}", shell=True)
    return pathf, (time() - startt)


def bch_dist_mash(pathf, threads, distdest=None):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    startt = time()
    subprocess.check_call(f"mash dist -p {threads} {pathf} > {distdest}", shell=True)
    return pathf, (time() - startt)


def bch_sketch_bindash(pathf, k, threads, *, destp, bbits, size):
    startt = time()
    subprocess.check_call(f"bindash sketch --kmerlen={k} --sketchsize64={size//64} --listfname={pathf} --outfname={destp} --nthreads={threads}", shell=True)
    return destp, (time() - startt)


def repeat_x(func, ntimes, *args, **kwargs):
    timevec = np.zeros(ntimes)
    for i in range(ntimes):
        res = func(*args, **kwargs)
        timevec[i] = res[-1]
    return (res[0], np.median(timevec))


def test():
    print(bch_sketch_mash("fn.txt", 31, 1, "mashdest.msh", 1024))
    bch_dist_mash("mashdest.msh", 1, "mashDistance.msh")
    bdpath, bdtime = bch_sketch_bindash("fn.txt", 31, 1, size=1024, destp="bdshdest.bdsh", bbits=32)
    bch_dist_bindash(bdpath, 1, 1024, "bdshDistance.bdsh")
    print(bch_dist_dashing2("fn.txt", 31, 1, 1024, distdest="d1dist.bin", binary=True, regsize=8))
    print(bch_dist_dashing("fn.txt", 31, 1, 1024, distdest="d2dist.bin", binary=True))
    print(bch_sketch_mash("fn.txt", 31, 1, "mashdest.msh", 1024))
    print(bch_sketch_dashing("fn.txt", 31, 1, 1024))
    print(bch_sketch_dashing2("fn.txt", 31, 1, 1024))
    print(bch_sketch_dashing2("fn.txt", 31, 1, 1024, oneperm=True))
    print(bdpath, bdtime)
    mmed = repeat_x(bch_sketch_mash, 3, "fn.txt", 31, 1, "mashdest.msh", 1024)
    print(mmed)
    return 0


def main():
    from argparse import ArgumentParser as AP
    ap = AP()
    ap.add_argument("--nrepeat", type=int, default=3)
    ap.add_argument("--nthreads", type=int, default=1)
    ap.add_argument("--sketchsize", type=int, default=0)
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("fnames", type=str)
    args = ap.parse_args()
    ssz = args.sketchsize
    if ssz <= 0:
        raise ValueError("sketchsize must be > 0")
    elif (ssz - 1) & ssz:
        raise ValueError("Sketchsize must be a power of 2 for this experiment")
    nt = args.nthreads
    fn = args.fnames
    assert os.path.isfile(fn), f"fnames argument ({fn}) does not exist"
    ntimes = args.nrepeat
    print(f"##Results for {fn}")
    print("#Method\tk\tSketchTime\tDistanceTime")
    mdfile = "mashdest.msh"
    mdstfile = "mashdist.out"
    msout_fn, tsketch = repeat_x(bch_sketch_mash, args.nrepeat, fn, args.k, nt, mdfile, size=args.sketchsize)
    msdistout_fn, tdist = repeat_x(bch_dist_mash, args.nrepeat, msout_fn, nt, mdstfile)
    print(f"Mash\t{args.k}\t{tsketch}\t{tdist}")
    return 0

if __name__ == "__main__":
    if not sys.argv[1:]:
        sys.exit(test())
    else:
        sys.exit(main())
