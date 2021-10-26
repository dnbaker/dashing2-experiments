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
    ap = AP("Parser taking a set of genomes, a set of ids, and makes a Makefile which generates the final table.", epilog="You will need to compile `safecat`, which cats and annotates the lines compiled.")
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
    subcombs = list(itertools.product((1<<x for x in sizeset), k))
    tupcombs = reduce(lambda x, y: x + y, ([(x, (y, z)) for y, z in subcombs] for x in tups))
    ruletups = [(f"computeonerow.py {full_genome_ids[lid]} {full_genome_ids[rid]} {sz} {k} > {lid}.{rid}.{sz}.{k}.{outpref}.row", f"{lid}.{rid}.{sz}.{k}.{outpref}.row", (full_genome_ids[lid], full_genome_ids[rid]), sz, k) for (lid, rid), (sz, k) in tupcombs]
    outf = f"outfiles.{outpref}.txt"
    with open(outf, "w") as f:
        for t in ruletups:
            print(t[1], file=f)
    with open(ap.outfile, "w") as f:
        f.write(f"all: experiment_result.{outpref}\n")
        f.write(f"experiment_result.{outpref}: {' '.join(x[1] for x in ruletups)}\n\tsafecat.py {outf} experiment_result.{outpref}\n")
        for line in itertools.starmap(makerule, (x[:2] for x in ruletups)):
            print(line, file=f)
