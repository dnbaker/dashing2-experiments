from sketchingfns import *

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
    import os
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
                # For now, only compute PMH
                for isbin, bstr in zip((True, False), ("-bin", "-txt")):
                    for regsize in (8, 4, 2, 1, .5):
                        for cssize in [500000, 2500000, 10000000]:
                            OP2 = "PMH" + bstr
                            OP3 = OP2 + "-%g-%s" % (regsize, str(cssize) if cssize is not None else "exact")
                            distdest = f"d2dest.{OP3}.k{k}.{rstr}"
                            d2out_fn, tsketch = repeat_x(bch_sketch_pmh, args.nrepeat, fn, k=k, threads=nt, size=ssz, cssize=cssize)
                            d2distout_fn, tdist = repeat_x(bch_dist_pmh, args.nrepeat, fn, k=k, threads=nt, size=ssz, cssize=cssize, regsize=regsize, binary=isbin, distdest=distdest)
                            print(f"{OP3}\t{k}\t{ssz}\t{regsize}\t{nt}\t{tsketch}\t{tdist}", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
