import sys

import numpy as np

import sketch
import multiprocessing as mp
from argparse import ArgumentParser as ArgP


ap = ArgP("Screen simulator")
aa = ap.add_argument
aa("--ngenomes", '-G', type=int, default=10000)
aa("--setsize", '-Z', type=int, default=2500000)
aa("--startsketchsize", '-s', type=int, default=8)
aa("--stopsketchsize", '-S', type=int, default=16)
aa("--tmpfileprefix", "-T", type=str, default="simpref")
aa("--samplefrac", "-F", type=float, default=0.1)
args = ap.parse_args()

prefix = args.tmpfileprefix

ssrange = range(args.startsketchsize, args.stopsketchsize)


setsizes = np.random.poisson(args.setsize, args.ngenomes)
sets = list(map(sketch.util.randset, setsizes))
nsampled = np.array([int(np.ceil(args.samplefrac * len(x))) for x in sets])
uset = np.hstack([x[:ns] for x, ns in zip(sets, nsampled)] + [sketch.util.randset(args.setsize * 5)])
total_items = sum(map(len, sets))
print("MeanAbsErr\tMeanAbsBias\tMaxErr\tMeanAbsRelErr\tMaxRelErr")
print(f"uset: {len(uset)}. All sets: {total_items}", file=sys.stderr)

for setsketchsize in range(args.startsketchsize, args.stopsketchsize):
    nregs = 1 << setsketchsize
    print("Processing for size = %d" % setsketchsize, file=sys.stderr)
    sketches = list(map(lambda x: sketch.setsketch.css_from_np(x, nregs), sets))
    mainsketch = sketch.setsketch.css_from_np(uset, nregs)
    usizes = np.array([len(uset) + len(oset) - ns for oset, ns in zip(sets, nsampled)])
    jaccards = nsampled / usizes
    print("Mean true jaccard: %g" % np.mean(jaccards))
    cjaccs = np.array([mainsketch.jaccard_index(rhs) for rhs in sketches])
    mabserr = np.mean(np.abs(cjaccs - jaccards))
    mabsrelerr = np.mean(np.abs((cjaccs - jaccards) / jaccards))
    mabsbias = np.mean(cjaccs - jaccards)
    maxerr = np.max(np.abs(cjaccs - jaccards))
    maxrelerr = np.max((np.abs(cjaccs - jaccards) / jaccards))
    print("%0.10g\t%0.10g\t%0.10g\t%0.10g\t%0.10g\n" % (mabserr, mabsbias, maxerr, mabsrelerr, maxrelerr))
