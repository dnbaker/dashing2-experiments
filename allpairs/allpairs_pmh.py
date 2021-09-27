import sys
import os
import numpy as np
import subprocess
from time import time


def bch_sketch_mash(pathf, k, threads, destp, size):
    startt = time()
    subprocess.check_call(f"mash sketch -k {k} -o {destp} -l {pathf} -p {threads}", shell=True)
    stopt = time()
    destp += ".msh"
    assert os.path.isfile(destp)
    return destp, (stopt - startt)


def bch_sketch_dashing(pathf, k, threads, size):
    startt = time()
    subprocess.check_call(f"dashing sketch -S {int(np.ceil(np.log2(size)))} -k {k} -F {pathf} -p {threads}, 2>/dev/null", shell=True)
    stopt = time()
    return pathf, (stopt - startt)


def bch_dist_dashing(pathf, k, threads, size, distdest=None, binary=False):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing")
    startt = time()
    pstr = " --emit-binary" if binary else ""
    subprocess.check_call(f"dashing dist --cache-sketches -k {k} -O{distdest} -o {distdest + '.sizes'} -S {int(np.ceil(np.log2(size)))} -F {pathf} -p {threads} {pstr} 2>/dev/null", shell=True)
    return pathf, (time() - startt)


def bch_sketch_dashingbb(pathf, k, threads, size):
    startt = time()
    subprocess.check_call(f"dashing sketch --use-bb-minhash --bbits 8 -S {int(np.ceil(np.log2(size)))} -k {k} -F {pathf} -p {threads}, 2>/dev/null", shell=True)
    stopt = time()
    return pathf, (stopt - startt)


def bch_dist_dashingbb(pathf, k, threads, size, distdest=None, binary=False):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashingbb")
    startt = time()
    pstr = " --emit-binary" if binary else ""
    subprocess.check_call(f"dashing dist --use-bb-minhash --bbits 8 --cache-sketches -k {k} -O{distdest} -o {distdest + '.sizes'} -S {int(np.ceil(np.log2(size)))} -F {pathf} -p {threads} {pstr} 2>/dev/null", shell=True)
    return pathf, (time() - startt)


def bch_sketch_pmh(pathf, k, threads, size, *, cssize=None):
    startt = time()
    csstr = f"--countsketch-size {cssize}" if cssize is not None else ""
    subprocess.check_call(f"dashing2 sketch --probminhash -k {k} -S {size} -F {pathf} -p {threads} {csstr} 2>/dev/null ", shell=True)
    return pathf, (time() - startt)


def bch_dist_pmh(pathf, k, threads, size, *, cssize=None, distdest=None, binary=False, regsize=-1):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_pmh")
    if regsize <= 0. or regsize not in [0.5, 1, 2, 4, 8]:
        raise ValueError("regsize must be > 0 and in [.5, 1, 2, 4, 8]")
    startt = time()
    pstr = " --binary-output " if binary else ""
    pstr += f" --countsketch-size {cssize}" if cssize is not None else ""
    subprocess.check_call(f"dashing2 sketch -k {k} --bbit-sigs --fastcmp {regsize} --cache --cmpout {distdest} -S {size} -F {pathf} -p {threads} {pstr} 2>/dev/null", shell=True)
    return pathf, (time() - startt)



def bch_sketch_dashing2(pathf, k, threads, size, oneperm=False):
    startt = time()
    pstr = " --oph " if oneperm else ""
    subprocess.check_call(f"dashing2 sketch -k {k} -S {size} -F {pathf} -p {threads} {pstr} 2>/dev/null ", shell=True)
    return pathf, (time() - startt)


def bch_dist_dashing2(pathf, k, threads, size, oneperm=False, distdest=None, binary=False, regsize=-1):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    if regsize <= 0. or regsize not in [0.5, 1, 2, 4, 8]:
        raise ValueError("regsize must be > 0 and in [.5, 1, 2, 4, 8]")
    startt = time()
    pstr = " --oph " if oneperm else ""
    pstr += " --binary-output " if binary else ""
    subprocess.check_call(f"dashing2 sketch -k {k} --fastcmp {regsize} --cache --cmpout {distdest} -S {size} -F {pathf} -p {threads} {pstr} 2>/dev/null", shell=True)
    return pathf, (time() - startt)


def bch_dist_bindash(pathf, threads, distdest=None):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    startt = time()
    subprocess.check_call(f"bindash dist --outfname={distdest} --nthreads={threads} {pathf} 2>/dev/null", shell=True)
    return pathf, (time() - startt)


def bch_dist_mash(pathf, threads, distdest=None):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    startt = time()
    subprocess.check_call(f"mash triangle -p {threads} {pathf} > {distdest}", shell=True)
    return pathf, (time() - startt)


def bch_sketch_bindash(pathf, k, threads, *, destp, bbits, size):
    startt = time()
    subprocess.check_call(f"bindash sketch --kmerlen={k} --sketchsize64={size//64} --listfname={pathf} --outfname={destp} --nthreads={threads} 2>/dev/null", shell=True)
    return destp, (time() - startt)


def repeat_x(func, ntimes, *args, **kwargs):
    timevec = np.zeros(ntimes)
    for i in range(ntimes):
        res = func(*args, **kwargs)
        timevec[i] = res[-1]
    return (res[0], np.median(timevec))


def test():
    print(bch_sketch_mash("fn.txt", 31, 1, "MASHdest", 1024))
    bch_dist_mash("MASHdest.msh", 1, "mashDistance.msh")
    bdpath, bdtime = bch_sketch_bindash("fn.txt", 31, 1, size=1024, destp="bdshdest.bdsh", bbits=32)
    bch_dist_bindash(bdpath, 1, "bdshDistance.bdsh")
    print(bch_dist_dashing2("fn.txt", 31, 1, 1024, distdest="d1dist.bin", binary=True, regsize=8))
    print(bch_dist_dashing("fn.txt", 31, 1, 1024, distdest="d2dist.bin", binary=True))
    print(bch_sketch_mash("fn.txt", 31, 1, "MASHdest", 1024))
    print(bch_sketch_dashing("fn.txt", 31, 1, 1024))
    print(bch_sketch_dashing2("fn.txt", 31, 1, 1024))
    print(bch_sketch_dashing2("fn.txt", 31, 1, 1024, oneperm=True))
    print(bch_sketch_pmh("fn.txt", 31, 1, 1024, cssize=100000))
    print(bch_dist_pmh("fn.txt", 31, 1, 1024, regsize=2, cssize=100000, distdest="pmhdest.bin", binary=True))
    print(bdpath, bdtime)
    mmed = repeat_x(bch_sketch_mash, 3, "fn.txt", 31, 1, "MASHdest", 1024)
    print(mmed)
    return 0


def main():
    from argparse import ArgumentParser as AP
    ap = AP()
    ap.add_argument("--nrepeat", type=int, default=3)
    ap.add_argument("--nthreads", '-p', action='append')
    ap.add_argument("--sketchsize", '-s', action='append')
    ap.add_argument("-k", action='append')
    ap.add_argument("fnames", type=str)
    args = ap.parse_args()
    if not args.nthreads:
        args.nthreads = [-1]
    fn = args.fnames
    print(f"##Options: {args}")
    print(f"##Results for {fn}")
    print("#Method\tk\tNumReg\tRegisterSize\tNumThreads\tSketchTime\tDistanceTime", flush=True)
    if not args.sketchsize:
        raise ValueError("Required: at least one value for sketchsize")
    if not args.k:
        raise ValueError("Required: at least one value for k")
    sszes = list(map(int, args.sketchsize))
    if any(x <= 0 or ((x - 1) & x) != 0 for x in sszes):
        raise ValueError("sketchsize must be > 1 and a power of 2")
    nt = args.nthreads
    assert os.path.isfile(fn), f"fnames argument ({fn}) does not exist"
    ntimes = args.nrepeat
    import random
    rstr = "".join(random.choice("abcdefgh") for i in range(12))
    for nt in map(int, args.nthreads):
        if nt < 0:
            from multiprocessing import cpu_count
            nt = cpu_count()
        for k in map(int, args.k):
            for ssz in sszes:
                # Handle MASH
                mdfile = f"MASHdest.k{k}.sz{ssz}.{rstr}"
                mdstfile = f"MASHdist.k{k}.sz{ssz}.{rstr}.phylip"
                msout_fn, tsketch = repeat_x(bch_sketch_mash, args.nrepeat, fn, k, nt, mdfile, size=ssz)
                msdistout_fn, tdist = repeat_x(bch_dist_mash, args.nrepeat, msout_fn, nt, mdstfile)
                print(f"Mash\t{k}\t{ssz}\t8\t{nt}\t{tsketch}\t{tdist}", flush=True)
                # Handle BinDash
                for regsize in (8, 4, 2, 1, .5):
                    nbits = int(regsize * 8)
                    bdfile = f"BDASHdest.k{k}.sz{ssz}.{rstr}.{nbits}"
                    bdstfile = f"BDASHdist.k{k}.sz{ssz}.{rstr}.{nbits}.out"
                    bdsout_fn, tsketch = repeat_x(bch_sketch_bindash, args.nrepeat, fn, k, threads=nt, destp=bdfile, bbits=nbits, size=ssz)
                    bddistout_fn, tdist = repeat_x(bch_dist_bindash, args.nrepeat, bdsout_fn, threads=nt, distdest=bdstfile)
                    print(f"Bindash-{regsize}\t{k}\t{ssz}\t{regsize}\t{nt}\t{tsketch}\t{tdist}", flush=True)
                for isbin in [1, 0]:
                    for isbb in [1, 0]:
                        sketchfn = bch_sketch_dashingbb if isbb else bch_sketch_dashing
                        distfn = bch_dist_dashingbb if isbb else bch_dist_dashing
                        d1out_fn, tsketch = repeat_x(sketchfn, args.nrepeat, fn, k=k, threads=nt, size=ssz)
                        d1distout_fn, tdist = repeat_x(distfn, args.nrepeat, fn, k=k, threads=nt, size=ssz, distdest=msdistout_fn.lower().replace("mash", "dashing"), binary=isbin)
                        basename = "Dashing1" if isbb else "Dashing1-bb"
                        print(f"{basename}-{'bin' if isbin else 'txt'}\t{k}\t{ssz}\t1\t{nt}\t{tsketch}\t{tdist}", flush=True)
                for isbin, bstr in zip((True, False), ("-bin", "-txt")):
                    for OP, D2S in zip([True, False], ("D2OP", "D2FSS")):
                        for regsize in (8, 4, 2, 1, .5):
                            OP2 = D2S + bstr
                            OP3 = OP2 + "-%g" % regsize
                            distdest = f"d2dest.{OP3}.k{k}.{rstr}"
                            d2out_fn, tsketch = repeat_x(bch_sketch_dashing2, args.nrepeat, fn, k=k, threads=nt, size=ssz, oneperm=OP)
                            d2distout_fn, tdist = repeat_x(bch_dist_dashing2, args.nrepeat, fn, k=k, threads=nt, size=ssz, oneperm=OP, regsize=regsize, binary=isbin, distdest=distdest)
                            print(f"{OP3}\t{k}\t{ssz}\t{regsize}\t{nt}\t{tsketch}\t{tdist}", flush=True)
                for isbin, bstr in zip((True, False), ("-bin", "-txt")):
                    for regsize in (8, 2, 1, .5):
                        for cssize in [500000, 2500000, None]:
                            OP2 = "PMH" + bstr
                            OP3 = OP2 + "-%g-%s" % (regsize, str(cssize) if cssize is not None else "exact")
                            distdest = f"d2dest.{OP3}.k{k}.{rstr}"
                            d2out_fn, tsketch = repeat_x(bch_sketch_pmh, args.nrepeat, fn, k=k, threads=nt, size=ssz, cssize=cssize)
                            d2distout_fn, tdist = repeat_x(bch_dist_pmh, args.nrepeat, fn, k=k, threads=nt, size=ssz, cssize=cssize, regsize=regsize, binary=isbin, distdest=distdest)
                            print(f"{OP3}\t{k}\t{ssz}\t{regsize}\t{nt}\t{tsketch}\t{tdist}", flush=True)
    return 0

if __name__ == "__main__":
    if not sys.argv[1:]:
        sys.exit(test())
    else:
        sys.exit(main())
