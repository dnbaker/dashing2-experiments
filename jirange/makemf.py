import sys
from functools import reduce
import string
import random
import itertools
from jir import getall, parsedata


def makerule(cmd, name):
    return f"{name}:\n\t{cmd}"


if __name__ == "__main__":
    from argparse import ArgumentParser as AP
    ap = AP()
    ap.add_argument("table", help="Path to a list of selected genomes with ranged Jaccards")
    ap.add_argument("fnames", help="Path to a list of selected genomes")
    ap.add_argument("--outfile", default="/dev/stdout", help="Path to use for output. Default: stdout")
    ap.add_argument("-s", type=int, default=10, help="Set start sketch size in log2.")
    ap.add_argument("-S", type=int, default=14, help="Set start sketch size in log2.")
    ap.add_argument("-T", type=int, default=1, help="Set step size for sketch size.")
    ap.add_argument("-k", type=int, action='append', help="Values of K to use.")
    ap = ap.parse_args()
    
    outpref = "".join(random.choice(string.ascii_letters) for x in range(10))
    if not ap.k:
        print("Using k = 31 by default", file=sys.stderr)
        ap.k = [31]
    k = ap.k
    
    res = parsedata(ap.table, ap.fnames)
    bkts, rngs, tups, total_ids, full_genome_ids = res
    tups = [tuple(map(int, x)) for x in tups.reshape(-1, 2)]
    selected = [full_genome_ids[i] for i in total_ids]
    tupids = [str(hash(full_genome_ids[i] + full_genome_ids[j])) for i, j in tups]
    sizeset = list(range(ap.s, ap.S, ap.T))
    typecombs = itertools.product((1<<x for x in sizeset), k)
    subcombs = list(itertools.product((1<<x for x in sizeset), k))
    tupcombs = reduce(lambda x, y: x + y, ([(x, (y, z)) for y, z in subcombs] for x in tups))
    ruletups = [(f"computeonerow.py {full_genome_ids[lid]} {full_genome_ids[rid]} {sz} {k} > {lid}.{rid}.{sz}.{k}.row", f"{lid}.{rid}.{sz}.{k}.row", (full_genome_ids[lid], full_genome_ids[rid]), sz, k) for (lid, rid), (sz, k) in tupcombs]
    rules = [x[1] for x in ruletups]
    with open("outfiles.txt", "w") as f:
        for r in rules:
            print(r, file=f)
    ruletexts = [makerule(cmd, dest) for cmd, dest, _, __, ___ in ruletups]
    with open(ap.outfile, "w") as f:
        f.write(f"all: experiment_result.{outpref}\n")
        f.write(f"experiment_result.{outpref}: {' '.join(rules)}\n\tsafecat.sh outfiles.txt experiment_result.{outpref}\n")
        for rt in ruletups:
            print(makerule(rt[0], rt[1]), file=f)
