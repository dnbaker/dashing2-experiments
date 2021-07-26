import numpy as np
import scipy.sparse as sp
import sys
from sys import stderr
import subprocess
import argparse
import os
from subprocess import PIPE, check_call


def check_output(x):
    from subprocess import check_output as sco
    t = 0
    while 1:
        try:
            return sco(x, shell=True, stderr=PIPE)
        except Exception as e:
            t += 1
            if t == 10:
                print("Failed to subprocess call '%s', error = %s" % (x, e), file=sys.stderr)
                raise


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
    print(f"Total {len(total_ids)} ids parsed", file=sys.stderr)
    return list(map(np.array, (bkts, rngs, tups, total_ids)))


def parsedata(path, fn):
    bkts, rngs, tups, ids = parse_bf(path)
    genome_ids = [x.strip() for x in open(fn)]
    selected_genomes = [genome_ids[i] for i in ids.ravel()]
    #print(len(selected_genomes))
    return (bkts, rngs, tups, ids, genome_ids)


def getani(l, r):
    '''For genomes l and r, which must exist as files, estimate ANI using fastANI
       returns a numpy array with [ANI, num, denom] as values.
    '''
    # return [ANI, num, denom]
    cmd = f"fastANI --minFraction 0. -q {l} -r {r} -o /dev/stdout"
    out = check_output(cmd).decode().strip()
    try:
        return np.fromiter(out.split("\t")[-3:], dtype=np.float64)[0]
    except Exception as se:
        return 0.


def getmashji(l, r, *, k, size=1024):
    return float(check_output(f"jaccard_mash dist -s {size} -k {k} -j -t {l} {r}").decode().strip().split("\n")[1].split()[-1])


def getdashingji(l, r, *, k, l2s=10):
    return float(check_output(f"dashing dist -S{l2s} -k{k} {l} {r}").decode().strip().split("\n")[-2].split("\t")[-1])


def exact_wjaccard(p1, p2, *, k):
    return float(check_output(f"dashing2 sketch -k{k} --countdict --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def exact_jaccard(p1, p2, *, k):
    res = check_output(f"dashing2 sketch -k{k} --set --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n')
    rf = float(res.split('\n')[-2].split("\t")[1])
    if rf > 1.:
        print("Res", res, rf, file=sys.stderr)
    return rf


def setsketch_jaccard(p1, p2, size, *, k, nb=8, fss=False, executable="dashing2"):
    '''k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       nb is 8 by default, which is Dashing 2's typical size, though dashing-f and dashing-ld use 4 bytes and 16 bytes
       to store float32 and long double hash registers, respectively.
    '''
    fcstr = f"--fastcmp {nb}" if nb < 8 else ""
    return float(check_output(f"{executable} sketch -k{k} {'--full-setsketch ' if fss else ''} -S{size} --fastcmp {nb} --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def bbminhash_jaccard(p1, p2, size, *, k, nb=32, fss=False, executable="dashing2"):
    '''k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       nb is 8 by default, which is Dashing 2's typical size
    '''
    assert (nb & (nb - 1)) == 0 and 4 <= nb <= 64, "nb, number of bits for bbit minhash, should be 4, 8, 16, 32, or 64"
    fcstr = f"--fastcmp {nb / 8}"
    return float(check_output(f"{executable} sketch -k{k} {'--full-setsketch ' if fss else ''} --bbit-sigs -S{size} {fcstr} --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def bindash_jaccard(p1, p2, size, *, k, nb=8, executable="bindash"):
    """Get bindash Jaccard, as well as caching its sketches
    """
    assert size % 64 == 0, "Size must be divisible by 64"
    bb = nb * 8
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
        cmd = f"{executable} dist {cp1} {cp2}"
        # print(f"bindash failed to emit a result {cmd}. Returning -1; this signals an error.", file=sys.stderr)
        return -1
    toks = out.split('\t')[-1].split("/")
    # print(toks)
    num, denom = map(float, toks)
    # print(num, denom)
    return num / denom


def bagminhash_jaccard(p1, p2, size, *, k, nb=8, cssize=-1):
    css = "--countsketch-size %d" % cssize if cssize > 0 else ""
    fss = f" --fastcmp {nb}" if nb in (8, 4, 2, 1, .5) else ""
    cstr = f"dashing2 sketch --cache --cmpout /dev/stdout -k {k} --multiset {css + fss} {p1} {p2}"
    return check_output(cstr).decode().strip('\n').split('\n')[-2].split("\t")[1]


def probminhash_jaccard(p1, p2, size, *, k, nb=8, cssize=-1):
    css = "--countsketch-size %d" % cssize if cssize > 0 else ""
    fss = f" --fastcmp {nb}" if nb in (8, 4, 2, 1, .5) else ""
    cstr = f"dashing2 sketch --cache --cmpout /dev/stdout -k {k} --prob {css + fss} {p1} {p2}"
    return check_output(cstr).decode().strip('\n').split('\n')[-2].split("\t")[1]


header = "#G1\tG2\tK\tsketchsize\tANI\tWJI\tJI\tMash\tDash1\tBD8\tBD4\tBD2\tBD1\tBDN\tSS8\tSS2\tSS1\tSSN\tFSS8\tFSS2\tFSS1\tFSSN\tMH8\tMH4\tMH2\tMH1\tMHN\tFMH8\tFMH4\tFMH2\tFMH1\tFMHN"
PMNBs = [8, 4, 2, 1, .5]
BMNBs = [8, 4, 2, 1, .5]
CSSZ = [-1, 50000000, 50000]
PMHSettings = [(b, cs) for b in PMNBs for cs in CSSZ]
BMHSettings = [(b, cs) for b in BMNBs for cs in CSSZ]
for (b, cs) in PMHSettings:
    header = header + "\tPMH%s%s" % (b if b >= 1 else "N", "-%d" % cs if cs > 0 else "Exact")
for (b, cs) in BMHSettings:
    header = header + "\tBMH%s%s" % (b if b >= 1 else "N", "-%d" % cs if cs > 0 else "Exact")


def getall(l, r, k=17, size=1024, executable="dashing2", cssize_set=[500, 50000, 50000000]):
    '''
        for values of k, size, and executable, return
        all similarity comparisons using Mash, Dashing, Dashing2, and fastANI
    '''
    ret = np.array([getani(l, r),
                     exact_wjaccard(l, r, k=k),
                     exact_jaccard(l, r, k=k),
                     getmashji(l, r, k=k, size=size),
                     getdashingji(l, r, k=k, l2s=int(np.log2(size)))] +
                     [bindash_jaccard(l, r, k=k, size=size, nb=nb) for nb in (8, 4, 2, 1, .5)] +
                     [setsketch_jaccard(l, r, size=size, k=k, nb=nb, fss=False, executable=executable) for nb in (8, 2, 1, .5)] +
                     [setsketch_jaccard(l, r, size=size, k=k, nb=nb, fss=True, executable=executable) for nb in (8, 2, 1, .5)] +
                     [bbminhash_jaccard(l, r, size=size, k=k, nb=int(nb * 8), fss=False, executable=executable) for nb in (8, 4, 2, 1, .5)] +
                     [bbminhash_jaccard(l, r, size=size, k=k, nb=int(nb * 8), fss=True, executable=executable) for nb in (8, 4, 2, 1, .5)] +
                     [probminhash_jaccard(l, r, size=size, k=k, nb=b, cssize=csz) for b, csz in PMHSettings] +
                     [bagminhash_jaccard(l, r, size=size, k=k, nb=b, cssize=csz) for b, csz in BMHSettings], np.float32)
    if ret[2] > 1.:
        print(f"Distance {l} {r}, k = {k}, size = {size}... For some reason, Jaccard was > 1...What's going on?", ret[2], file=sys.stderr)
    return ret

def packed(x):
    l, r, k, size, executable = x
    return getall(l, r, k=k, size=size, executable=executable)

def packedani(x):
    l, r = x
    return getani(l, r)


def pargetall(tups, k=17, executable="dashing2", cpu=-1):
    import multiprocessing as mp
    if cpu < 0: cpu = mp.cpu_count()
    with mp.Pool(cpu) as p:
        ret = np.stack(p.map(packed, tups))
        # print(ret)
        return ret



if __name__ == "__main__":
    if(sys.argv[2:] or "-h" in sys.argv or "--help" in sys.argv):
        ap = argparse.ArgumentParser()
        ap.add_argument("table", help="Path to a list of selected genomes with ranged Jaccards")
        ap.add_argument("fnames", help="Path to a list of selected genomes")
        ap.add_argument("-k", type=int, default=17, help="Set k for experiment. If exact representations aren't cached for this value of k, it may take a very long time to run")
        ap.add_argument("-s", type=int, default=10, help="Set start sketch size in log2.")
        ap.add_argument("-S", type=int, default=14, help="Set start sketch size in log2.")
        ap.add_argument("-T", type=int, default=2, help="Set step size for sketch size.")
        ap.add_argument("--cpu", type=int, default=-1)
        ap.add_argument("--executable", '-E', default="dashing2")
        ap.add_argument("--name", default="noname")
        ap.add_argument("-o", "--outfile", default="/dev/stdout")
        args = ap.parse_args()
        k = args.k
        res = parsedata(args.table, args.fnames)
        bkts, rngs, tups, total_ids, full_genome_ids = res
        selected = [full_genome_ids[i] for i in total_ids]
        print(len(selected), ", frac = %g" % (len(selected) / len(full_genome_ids)), file=sys.stderr)
        s = ""
        for tup in tups:
            # print(tup.shape)
            s += ",".join("%s-%s" % (x[0], x[1]) for x in tup)
            s += '\n'
        hv = abs(hash(",".join(sys.argv)))
        rng = range(args.s, args.S, args.T)
        sdict = {"k": k, "executable": args.executable, "cpu": args.cpu}
        tdict = {i: (full_genome_ids[l], full_genome_ids[r]) for i, (l, r) in enumerate(tups.reshape(-1, 2))}
        tups = [(full_genome_ids[l], full_genome_ids[r], k, 1 << size, args.executable) for l, r in tups.reshape(-1, 2) for size in rng]
        with open("settings.%s.%d.txt" % (args.name, hv), "w") as f:
            for st in tups:
                l, r, k, size, _ = st
                print("%s\t%s\t%d\t%d" % (l, r, k, size), file=f)
        print(f"Generated {len(tups)} tuples, which are being passed to pargetall", file=sys.stderr)
        fullmat = np.stack(pargetall(tups, **sdict))
        print(fullmat.shape)
        fs = str(fullmat.shape).replace(" ", "").replace(",", "-")
        fullmat.astype(np.float32).tofile("fullmat.%s.f32.%d.%s" % (args.name, hv, fs))
        with open(args.outfile, "w") as ofp:
            print(header, file=ofp)
            for (l, r, k, size, _), mr in zip(tups, fullmat):
                    print(f"{l}\t{r}\t{k}\t{size}\t" + "\t".join(map(str, mr)), file=ofp)
    else:
        print("Running tests, not running experiment", file=sys.stderr)
        parse_bf(sys.argv[1] if sys.argv[1:] else "selected_buckets_100.txt")
