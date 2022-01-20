from sketchingfns import *
import sys
import os
import glob


def test():
    print(bch_sketch_mash("fn.txt", 31, 1, "MASHdest.msh", 1024))
    bch_dist_mash("MASHdest.msh", 1, "mashDistance.msh")
    bdpath, bdtime = bch_sketch_bindash("fn.txt", 31, 1, size=1024, destp="bdshdest.bdsh", bbits=32)
    bch_dist_bindash(bdpath, 1, "bdshDistance.bdsh")
    print(bch_dist_dashing2("fn.txt", 31, 1, 1024, distdest="d1dist.bin", binary=True, regsize=8))
    print(bch_dist_dashing("fn.txt", 31, 1, 1024, distdest="d2dist.bin", binary=True))
    print(bch_sketch_mash("fn.txt", 31, 1, "MASHdest.msh", 1024))
    print(bch_sketch_dashing("fn.txt", 31, 1, 1024))
    print(bch_sketch_dashing2("fn.txt", 31, 1, 1024))
    print(bch_sketch_dashing2("fn.txt", 31, 1, 1024, oneperm=True))
    print(bdpath, bdtime)
    mmed = repeat_x(bch_sketch_mash, 3, "fn.txt", 31, 1, "MASHdest.msh", 1024)
    print(mmed)
    return 0

acceptable_types = ("mash", "bindash", "d2op", "d2fss", "d1", "pmh", None)


def fallible_rm(x):
    return

def main():
    from argparse import ArgumentParser as AP
    ap = AP()
    ap.add_argument("--nrepeat", type=int, default=3)
    ap.add_argument("--nthreads", '-p', action='append', type=int)
    ap.add_argument("--sketchsize", '-s', action='append', type=int)
    ap.add_argument("-k", action='append', type=int)
    ap.add_argument("--only-type", choices=acceptable_types, help="If set, only performs benchmark for a particular type.")
    ap.add_argument("--only-cssize", type=int, help="Specify this argument to override the count-sketch sizes for PMH sketching.", action='append')
    ap.add_argument("--pmhbbitsigs", action='store_true', help="Use this flag to use bbit sigs instead of processing as weighted SetSketch registers.")
    ap.add_argument("--executable", "--d2", default="dashing2", help="Alternate path to Dashing2 for varying architectures.")
    ap.add_argument("fnames", type=str)
    args = ap.parse_args()
    if not args.nthreads:
        args.nthreads = [-1]
    fn = args.fnames
    print(f"##Options: {args}")
    print(f"##Results for {fn}")
    print("#Method\tk\tNumReg\tRegisterSize\tNumThreads\tSketchTime\tDistanceTime", flush=True)
    def passes(key):
        assert key in acceptable_types
        return args.only_type in (key, None)

    if not args.sketchsize:
        raise ValueError("Required: at least one value for sketchsize")
    if not args.k:
        raise ValueError("Required: at least one value for k")
    sszes = args.sketchsize
    if any(x <= 0 or ((x - 1) & x) != 0 for x in sszes):
        raise ValueError("sketchsize must be > 1 and a power of 2")
    assert os.path.isfile(fn), f"fnames argument ({fn}) does not exist"
    ntimes = args.nrepeat
    import random
    rstr = "".join(random.choice("abcdefghijkl") for i in range(12))
    for nt in args.nthreads:
        if nt < 0:
            from multiprocessing import cpu_count
            nt = cpu_count()
        for k in args.k:
            for ssz in sszes:
                # This ternary skips dashing if not specified
                binvals = (1, 0) if passes("d1") else ()
                for isbin in binvals:
                    d1out_fn, tsketch = repeat_x(bch_sketch_dashing, args.nrepeat, fn, k=k, threads=nt, size=ssz)
                    d1distout_fn, tdist = repeat_x(bch_dist_dashing, args.nrepeat, fn, k=k, threads=nt, size=ssz, distdest=f"dashing1.dist.k{k}.sz{ssz}.{rstr}.phylip", binary=isbin)
                    print(f"Dashing1-{'bin' if isbin else 'txt'}\t{k}\t{ssz}\t1\t{nt}\t{tsketch}\t{tdist}", flush=True)
                    fallible_rm(d1distout_fn)
                if passes("d2op"):
                    for isbin, bstr in zip((True, False), ("-bin", "-txt")):
                        OP = True
                        D2S = "D2OP"
                        for regsize in (8, 4, 2, 1, .5):
                            OP2 = D2S + bstr
                            OP3 = OP2 + "-%g" % regsize
                            distdest = f"d2dest.{OP3}.k{k}.{rstr}"
                            d2out_fn, tsketch = repeat_x(bch_sketch_dashing2, args.nrepeat, fn, k=k, threads=nt, size=ssz, oneperm=OP, executable=args.executable)
                            d2distout_fn, tdist = repeat_x(bch_dist_dashing2, args.nrepeat, fn, k=k, threads=nt, size=ssz, oneperm=OP, regsize=regsize, binary=isbin, distdest=distdest, executable=args.executable)
                            print(f"{OP3}\t{k}\t{ssz}\t{regsize}\t{nt}\t{tsketch}\t{tdist}", flush=True)
                if passes("d2fss"):
                    for isbin, bstr in zip((True, False), ("-bin", "-txt")):
                        OP = False
                        D2S = "D2FSS"
                        for regsize in (8, 4, 2, 1, .5):
                            OP2 = D2S + bstr
                            OP3 = OP2 + "-%g" % regsize
                            distdest = f"d2dest.{OP3}.k{k}.{rstr}"
                            d2out_fn, tsketch = repeat_x(bch_sketch_dashing2, args.nrepeat, fn, k=k, threads=nt, size=ssz, oneperm=OP, executable=args.executable)
                            d2distout_fn, tdist = repeat_x(bch_dist_dashing2, args.nrepeat, fn, k=k, threads=nt, size=ssz, oneperm=OP, regsize=regsize, binary=isbin, distdest=distdest, executable=args.executable)
                            print(f"{OP3}\t{k}\t{ssz}\t{regsize}\t{nt}\t{tsketch}\t{tdist}", flush=True)
                            fallible_rm(d2distout_fn)
                # Handle MASH
                if passes("mash"):
                    mdfile = f"MASHdest.k{k}.sz{ssz}.{rstr}"
                    mdstfile = f"MASHdist.k{k}.sz{ssz}.{rstr}.phylip"
                    msout_fn, tsketch = repeat_x(bch_sketch_mash, args.nrepeat, fn, k, threads=nt, destp=mdfile, size=ssz)
                    msdistout_fn, tdist = repeat_x(bch_dist_mash, args.nrepeat, msout_fn, threads=nt, distdest=mdstfile)
                    print(f"Mash\t{k}\t{ssz}\t8\t{nt}\t{tsketch}\t{tdist}", flush=True)
                    set(map(fallible_rm, (msout_fn, msdistout_fn)))
                # Handle BinDash
                # This ternary skips bindash if not specified
                bdregsizes = (8, 4, 2, 1, .5) if passes("bindash") else ()
                for regsize in bdregsizes:
                    nbits = int(regsize * 8)
                    bdfile = f"BDASHdest.k{k}.sz{ssz}.{rstr}.{nbits}"
                    bdstfile = f"BDASHdist.k{k}.sz{ssz}.{rstr}.{nbits}.out"
                    bdsout_fn, tsketch = repeat_x(bch_sketch_bindash, args.nrepeat, fn, k, threads=nt, destp=bdfile, bbits=nbits, size=ssz)
                    bddistout_fn, tdist = repeat_x(bch_dist_bindash, args.nrepeat, bdsout_fn, threads=nt, distdest=bdstfile)
                    print(f"Bindash-{regsize}\t{k}\t{ssz}\t{regsize}\t{nt}\t{tsketch}\t{tdist}", flush=True)
                    fallible_rm(bddistout_fn)
                    set(map(fallible_rm, glob.iglob(f"{bdfile}*")))
                if passes("pmh"):
                    cssizes = (500000, 2500000) if args.only_cssize is None else tuple(args.only_cssize)
                    for isbin, bstr in zip((True, False), ("-bin", "-txt")):
                        for regsize in (8, 4, 2, 1, .5):
                            for cssize in cssizes:
                                OP2 = "PMH" + bstr
                                OP3 = OP2 + "-%g-%s" % (regsize, str(cssize) if cssize is not None else "exact")
                                distdest = f"d2dest.{OP3}.k{k}.{rstr}"
                                d2out_fn, tsketch = repeat_x(bch_sketch_pmh, args.nrepeat, fn, k=k, threads=nt, size=ssz, cssize=cssize, executable=args.executable)
                                d2distout_fn, tdist = repeat_x(bch_dist_pmh, args.nrepeat, fn, k=k, threads=nt, size=ssz, cssize=cssize, regsize=regsize, binary=isbin, distdest=distdest, bbit_sigs=args.pmhbbitsigs, executable=args.executable)
                                print(f"{OP3}\t{k}\t{ssz}\t{regsize}\t{nt}\t{tsketch}\t{tdist}", flush=True)
                                fallible_rm(d2distout_fn)
    return 0

if __name__ == "__main__":
    sys.exit(main())
