import numpy as np
import scipy.sparse as sp
import sys
from sys import stderr
import subprocess
import argparse
import os
from subprocess import PIPE, check_call


def check_output(x):
    from subprocess import Popen, PIPE
    t = 0
    while 1:
        process = Popen(x, shell=True, stdout=PIPE, stderr=PIPE, executable="/bin/bash")
        out, err = process.communicate()
        if process.returncode:
            t += 1
            print("Failed call time #%d '%s', error = %s" % (t, x, err.decode()), file=sys.stderr)
            if t == 10:
                raise RuntimeError("Failed 10 times to run command " + x)
        else:
            return out


def parse_bf(path):
    import itertools
    cfi = itertools.chain.from_iterable
    bkts = []
    rngs = []
    tups = []
    for line in open(path):
        toks = line.strip().split()
        bkts.append(int(toks[1]))
        rngs.append(list(map(float, toks[2].strip("[]").split("->"))))
        tups.append([list(map(int, x.split(":"))) for x in toks[3:]])
    bkts = np.array(bkts)
    ml = max(map(len, tups))
    total_ids = np.array(sorted(set(cfi(cfi(tups)))))
    for i in range(len(tups)):
        if len(tups[i]) < ml:
            diff = ml - len(tups[i])
            tups[i] += [[0,0]] * diff
        assert len(tups[i]) == ml
        if len(rngs[i]) < 2:
            rngs[i].append(rngs[i][0])
    assert len(set(map(len, rngs))) == 1
    rngs = np.array(rngs)
    # print(f"Total {len(total_ids)} ids parsed", file=sys.stderr)
    return list(map(np.array, (bkts, rngs, tups, total_ids)))


def parsedata(path, fn):
    bkts, rngs, tups, ids = parse_bf(path)
    genome_ids = [x.strip() for x in open(fn)]
    selected_genomes = [genome_ids[i] for i in ids.ravel()]
    #print(len(selected_genomes))
    return (bkts, rngs, tups, ids, genome_ids)


def getani(left, r, ex="fastANI"):
    '''For genomes left and r, which must exist as files, estimate ANI using
       fastANI returns a numpy array with [ANI, num, denom] as values.
    '''
    # return [ANI, num, denom]
    cmd = f"{ex} --minFraction 0. -q {left} -r {r} -o /dev/stdout"
    out = check_output(cmd).decode().strip()
    try:
        return np.fromiter(out.split("\t")[-3:], dtype=np.float64)[0]
    except Exception as se:
        return 0.


def getmashji(left, r, *, k, size=1024):
    if k > 32:
        return 0. # Mash doesn't support long kmers
    out = check_output(f"mash dist -s {size} -k {k} {left} {r}").decode().strip().split("\n")[-1].split()[-1]
    num, denom = map(float, out.split("/"))
    if not denom:
        return 0.
    return num / denom


def getdashingji(left, r, *, k, l2s=10):
    cyclic_flag = "" if k <= 32 else " --use-cyclic-hash "
    cmd = f"dashing dist {cyclic_flag} -S{l2s} -k{k} {left} {r}"
    return float(check_output(cmd).decode().strip().split("\n")[-2].split("\t")[-1])


def exact_wjaccard(p1, p2, *, k):
    return float(check_output(f"dashing2 sketch --phylip -k{k} --countdict --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def exact_jaccard(p1, p2, *, k):
    res = check_output(f"dashing2 sketch --phylip -k{k} --set --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n')
    rf = float(res.split('\n')[-2].split("\t")[1])
    if rf > 1.:
        print("Res", res, rf, file=sys.stderr)
    return rf


def fsarg2str(x):
    if x:
        return ""
    else:
        return " --oneperm "


def setsketch_jaccard(p1, p2, size, *, k, nb=8, fss=False, executable="dashing2"):
    '''k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       nb is 8 by default, which is Dashing 2's typical size, though dashing-f and dashing-ld use 4 bytes and 16 bytes
       to store float32 and long double hash registers, respectively.
    '''
    fcstr = f"--fastcmp {nb}" if nb < 8 else ""
    return float(check_output(f"{executable} sketch --phylip -k{k} {fsarg2str(fss)} -S{size} --fastcmp {nb} --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def bbminhash_jaccard(p1, p2, size, *, k, nb=32, fss=False, executable="dashing2"):
    '''k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       nb is 8 by default, which is Dashing 2's typical size
    '''
    assert (nb & (nb - 1)) == 0 and 4 <= nb <= 64, "nb, number of bits for bbit minhash, should be 4, 8, 16, 32, or 64"
    fcstr = f"--fastcmp {nb / 8}"
    return float(check_output(f"{executable} sketch --phylip -k{k} {fsarg2str(fss)} --bbit-sigs -S{size} {fcstr} --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def bindash_jaccard(p1, p2, size, *, k, nb=8, executable="bindash"):
    """Get bindash Jaccard, as well as caching its sketches
    """
    assert size % 64 == 0, f"Size must be divisible by 64, found {size}"
    bb = nb * 8
    if not os.path.isdir("bdsh_tmp"):
        os.mkdir("bdsh_tmp")
    def mc(p):
        return "bdsh_tmp/" + os.path.basename(p + f".{k}.{bb}.{size}")
    cp1, cp2 = map(mc, (p1, p2))
    for cp, p in zip((cp1, cp2), (p1, p2)):
        if not (os.path.isfile(cp + ".dat") and os.path.isfile(cp + ".txt")):
            cmd = f"{executable} sketch --kmerlen={k} --bbits={bb} --minhashtype=2 --sketchsize64={size // 64} --outfname={cp} {p}"
            check_call(cmd, shell=True, stderr=PIPE)
    out = check_output(f"{executable} dist {cp1} {cp2}").decode().strip()
    # print("Output: ", out)
    if not out:
        cmd = f"{executable} dist --ithres=1 {cp1} {cp2}"
        return -1
    toks = out.split('\t')[-1].split("/")
    num, denom = map(float, toks)
    return num / denom


def bagminhash_jaccard(p1, p2, size, *, k, nb=8, cssize=-1):
    css = "--countsketch-size %d" % cssize if cssize > 0 else ""
    fss = f" --fastcmp {nb}" if nb in (8, 4, 2, 1, .5) else ""
    cstr = f"dashing2 sketch --phylip --cache --bbit-sigs --cmpout /dev/stdout -S {size} -k {k} --multiset {css + fss} {p1} {p2}"
    return check_output(cstr).decode().strip('\n').split('\n')[-2].split("\t")[1]


def probminhash_jaccard(p1, p2, size, *, k, nb=8, cssize=-1):
    css = "--countsketch-size %d" % cssize if cssize > 0 else ""
    fss = f" --fastcmp {nb}" if nb in (8, 4, 2, 1, .5) else ""
    cstr = f"dashing2 sketch --phylip --cache --bbit-sigs --cmpout /dev/stdout -S {size} -k {k} --prob {css + fss} {p1} {p2}"
    return check_output(cstr).decode().strip('\n').split('\n')[-2].split("\t")[1]


header = "#G1\tG2\tK\tsketchsize\tANI\tWJI\tJI\tMash\tDash1\tBD8\tBD4\tBD2\tBD1\tBDN\tSS8\tSS4\tSS2\tSS1\tSSN\tFSS8\tFSS4\tFSS2\tFSS1\tFSSN\tMH8\tMH4\tMH2\tMH1\tMHN\tFMH8\tFMH4\tFMH2\tFMH1\tFMHN"
PMNBs = [8, 4, 2, 1, .5]
BMNBs = [8, 4, 2, 1, .5]
CSSZ = [-1, 50000000, 500000]
PMHSettings = [(b, cs) for b in PMNBs for cs in CSSZ]
BMHSettings = [(b, cs) for b in BMNBs for cs in CSSZ]
for (b, cs) in PMHSettings:
    header = header + "\tPMH%s%s" % (b if b >= 1 else "N", "-%d" % cs if cs > 0 else "Exact")
for (b, cs) in BMHSettings:
    header = header + "\tBMH%s%s" % (b if b >= 1 else "N", "-%d" % cs if cs > 0 else "Exact")


def getall(l, r, k=17, size=1024, executable="dashing2", cssize_set=[500000, 5000000, 50000000], faex="fastANI"):
    '''
        for values of k, size, and executable, return
        all similarity comparisons using Mash, Dashing, Dashing2, and fastANI
    '''
    bbnbs = (8, 4, 2, 1, .5)
    ret = np.array([getani(l, r, ex=faex),
                    exact_wjaccard(l, r, k=k),
                    exact_jaccard(l, r, k=k),
                    getmashji(l, r, k=k, size=size),
                    getdashingji(l, r, k=k, l2s=int(np.log2(size)))] +
                    [bindash_jaccard(l, r, k=k, size=size, nb=nb) for nb in bbnbs] +
                    [setsketch_jaccard(l, r, size=size, k=k, nb=nb, fss=False, executable=executable) for nb in bbnbs] +
                    [setsketch_jaccard(l, r, size=size, k=k, nb=nb, fss=True, executable=executable) for nb in bbnbs] +
                    [bbminhash_jaccard(l, r, size=size, k=k, nb=int(nb * 8), fss=False, executable=executable) for nb in bbnbs] +
                    [bbminhash_jaccard(l, r, size=size, k=k, nb=int(nb * 8), fss=True, executable=executable) for nb in bbnbs] +
                    [probminhash_jaccard(l, r, size=size, k=k, nb=b, cssize=csz) for b, csz in PMHSettings] +
                    [bagminhash_jaccard(l, r, size=size, k=k, nb=b, cssize=csz) for b, csz in BMHSettings], np.float32)
    if ret[2] > 1.:
        print(f"Distance {l} {r}, k = {k}, size = {size}... For some reason, Jaccard was > 1...What's going on?", ret[2], file=sys.stderr)
    return ret


def packed(x):
    from time import time
    startt = time()
    try:
        l, r, k, size, executable = x
    except:
        print(len(x), file=sys.stderr)
        raise
    ret = getall(l, r, k=k, size=size, executable=executable)
    return ret


def packedani(x):
    l, r = x
    return getani(l, r)


if __name__ == "__main__":
    if(sys.argv[2:] or "-h" in sys.argv or "--help" in sys.argv):
        ap = argparse.ArgumentParser()
        ap.add_argument("table", help="Path to a list of selected genomes with ranged Jaccards")
        ap.add_argument("fnames", help="Path to a list of selected genomes")
        ap.add_argument("-k", type=int, default=17, help="Set k for experiment. If exact representations aren't cached for this value of k, it may take a very long time to run")
        ap.add_argument("-s", type=int, default=10, help="Set start sketch size in log2.")
        ap.add_argument("-S", type=int, default=14, help="Set start sketch size in log2.")
        ap.add_argument("-T", type=int, default=1, help="Set step size for sketch size.")
        ap.add_argument("--cpu", type=int, default=-1)
        ap.add_argument("--executable", '-E', default="dashing2")
        ap.add_argument("--name", default="noname")
        ap.add_argument("-o", "--outfile", default=sys.stdout)
        args = ap.parse_args()
        k = args.k
        res = parsedata(args.table, args.fnames)
        bkts, rngs, tups, total_ids, full_genome_ids = res
        selected = [full_genome_ids[i] for i in total_ids]
        print(len(selected), ", frac = %g" % (len(selected) / len(full_genome_ids)), file=sys.stderr)
        s = ""
        for tup in tups:
            s += ",".join("%s-%s" % (x[0], x[1]) for x in tup)
            s += '\n'
        hv = abs(hash(",".join(sys.argv)))
        rng = range(args.s, args.S, args.T)
        sdict = {"k": k, "executable": args.executable, "cpu": args.cpu}
        tdict = {i: (full_genome_ids[left], full_genome_ids[r]) for i, (left, r) in enumerate(tups.reshape(-1, 2))}
        tups = [(full_genome_ids[left], full_genome_ids[r], k, 1 << size, args.executable) for left, r in tups.reshape(-1, 2) for size in rng]
        with open("settings.%s.%d.txt" % (args.name, hv), "w") as f:
            for st in tups:
                left, r, k, size, _ = st
                print("%s\t%s\t%d\t%d" % (left, r, k, size), file=f)
        print(f"Generated {len(tups)} tuples", file=sys.stderr)
        import multiprocessing as mp
        cpu = args.cpu
        if cpu < 0:
            cpu = mp.cpu_count()
        shape = (len(tups), 62)
        fs = str(shape).replace(" ", "").replace(",", "-").replace("(", "_").replace(")", "_")
        nblocks = (len(tups) + 1023) >> 10
        CS = 1024
        subtups = [tups[x:x + CS] for x in map(lambda x: x * CS, range(nblocks))]
        ofp = open(args.outfile, "w") if args.outfile is not sys.stdout else sys.stdout
        print(header, file=ofp, flush=True)
        rawmat = []
        with mp.Pool(cpu) as p:
            for i, st in enumerate(subtups):
                from time import time
                startt = time()
                print("Started subgroup %d/%d" % (i, len(subtups)), file=sys.stderr, flush=True)
                res = np.stack(list(p.map(packed, st)))
                rawmat.append(res)
                for (left, r, k, size, _), mr in zip(st, res.reshape(-1, 60)):
                    print(f"{left}\t{r}\t{k}\t{size}\t" + "\t".join(map(str, mr)), file=ofp, flush=True)
        if ofp is not sys.stdout:
            ofp.close()
        np.vstack(rawmat).tofile(f"fullmat.{args.name}.f32.{hv}.{k}.{fs}")
    else:
        print("Running tests, not running experiment", file=sys.stderr)
        parse_bf(sys.argv[1] if sys.argv[1:] else "selected_buckets_100.txt")
