from multiprocessing import cpu_count
import numpy as np
import os, sys
import gzip
from timeit import default_timer as time


def getknn(args, *, exact):
    from subprocess import check_call
    ts = time()
    knntype = "exact" if exact else "LSH-assisted"
    outfile = f"{args.inputfile}.{knntype}.{args.alphabet}.{args.sketchsize}.time{ts}.topk.matrix"
    sketchfile = f"{args.inputfile}.{knntype}.{args.alphabet}.{args.sketchsize}.time{ts}.sketches.bin"
    alphstr = {"DNA": ""
             , "PROTEIN20": "--protein20"
             , "PROTEIN": "--protein"
             , "PROTEIN14": "--protein14"
             , "PROTEIN6": "--protein6"}[args.alphabet]
    cmd = f"{args.executable} sketch --cache -k {args.kmer_length} -o {sketchfile} --cmpout {outfile} -p{args.nthreads} -S{args.sketchsize} {alphstr} --nLSH {args.nLSH} --topk {args.nneighbors} --sig-ram-limit {args.ramlimit}"
    if args.mode == "parse-by-seq":
        cmd += " --parse-by-seq " + args.inputfile
    else:
        cmd += " -F " + args.inputfile
    if exact:
        cmd = "EXACT_KNN=1 " + cmd
    print("cmd: ", cmd)
    check_call(cmd, shell=True)
    t = time() - ts
    print(f"Took {t}s to perform {'exhaustive' if exact else 'LSH-assisted'} comparisons.")
    return outfile, t


def parse_knn(path):
    ifp = open(path)
    names = [x.split()[0] for x in ifp][1:]
    name2id = {name: idx for idx, name in enumerate(names)}
    ifp.seek(os.SEEK_SET)
    next(ifp)
    return [np.fromiter(map(name2id.__getitem__, (x.split(":")[0] for x in line.split()[1:])), np.uint32) for line in ifp]


def firstline(x):
    try:
        return next(gzip.open(x, "r"))
    except gzip.BadGzipFile:
        return next(open(x))


def getargs():
    from argparse import ArgumentParser as AP, RawDescriptionHelpFormatter
    ap = AP("KNN experiment", formatter_class=RawDescriptionHelpFormatter)
    aa = ap.add_argument
    aa("inputfile", help="Path to input file.")
    aa("--kmer-length", "-k", type=int, default=-1)
    aa("--alphabet", choices=("DNA", "PROTEIN", "PROTEIN14", "PROTEIN6"), default="DNA")
    aa("--mode", choices=("parse-by-seq", "fileset", None), default=None, help="Mode for experiment - parse-by-seq parses individual sequences from inputfile, or fileset parses sequence sets from files at each line. If unset, inferred.")
    aa("--nthreads", "-p", type=int, default=-1, help="Number of threads to usee. if <= 0, uses all available threads (%d)" % cpu_count())
    aa("--sketchsize", "-S", type=int, default=512)
    aa("--nneighbors", "--neighbors", "--n_neighbors", "-N", type=int, default=25)
    aa("--nLSH", type=int, default=2)
    aa("--executable", default="dashing2")
    aa("--ramlimit", "-R", type=int, default=20 << 30, help="Maximum memory to use for maintaining signature matrix. Lower this to move computation to disk.")

    args = ap.parse_args()
    DKL = {"DNA": 31, "PROTEIN": 14, "PROTEIN20": 14, "PROTEIN14": 16, "PROTEIN6": 24}
    args.kmer_length = args.kmer_length if args.kmer_length > 0 else DKL[args.alphabet]
    args.nthreads = args.nthreads if args.nthreads > 0 else cpu_count()
    if args.mode is None:
        args.mode = "fileset" if os.path.isfile(firstline(args.inputfile).strip()) else "parse-by-seq"
    return args

if __name__ == "__main__":
    args = getargs()
    lshmax, lshtime = getknn(args, exact=False)
    exactmat, exacttime = getknn(args, exact=True)
