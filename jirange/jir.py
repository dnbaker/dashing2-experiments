import numpy as np
import scipy.sparse as sp
import sys
import subprocess
import argparse


def check_output(x):
    from subprocess import check_output as sco
    return sco(x, shell=True)


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
    return np.fromiter(check_output(cmd).decode().strip().split("\t")[-3:], dtype=np.float64)


def getmashji(l, r, *, k, sz=1024):
    return float(check_output(f"jaccard_mash dist -s {sz} -k {k} -j -t {l} {r}").decode().strip().split("\n")[1].split()[-1])


def getdashingji(l, r, *, k, l2s=10):
    return float(check_output(f"dashing dist -S{l2s} -k{k}{l} {r}").decode().strip().split("\n")[-2].split("\t")[-1])


def exact_wjaccard(p1, p2, k=17):
    return float(check_output(f"dashing2 sketch -k{k} --countdict --exact-kmer-dist --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def exact_jaccard(p1, p2, k=17):
    return float(check_output(f"dashing2 sketch -k{k} --set --exact-kmer-dist --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def setsketch_jaccard(p1, p2, size, k=17, nb=8, fss=False, executable="dashing2"):
    '''k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       nb is 8 by default, which is Dashing 2's typical size, though dashing-f and dashing-ld use 4 bytes and 16 bytes
       to store float32 and long double hash registers, respectively.
    '''
    return float(check_output(f"{executable} sketch -k{k} {'--full-setsketch ' if fss else ''} -S{size} --fastcmp {nb} --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def bbminhash_jaccard(p1, p2, size, k=17, bbnb=32, fss=False, executable="dashing2"):
    '''k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       bbnb is 8 by default, which is Dashing 2's typical size
    '''
    nb = bbnb // 8
    assert (bbnb & (bbnb - 1)) == 0 and 4 <= bbnb <= 64, "bbnb, number of bits for bbit minhash, should be 4, 8, 16, 32, or 64"
    return float(check_output(f"{executable} sketch -k{k} {'--full-setsketch ' if fss else ''} --bbit-sigs -S{size} --fastcmp {nb} --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def getall(l, r, k=17, sz=1024, executable="dashing2"):
    '''
        for values of k, sz, and executable, return
        all similarity comparisons using Mash, Dashing, Dashing2, and fastANI
    '''
    return np.array([getmashji(l, r, k=k, sz=sz),
                     getdashingji(l, r, k=k, l2s=int(np.log2(sz))),
                     exact_jaccard(l, r),
                     exact_wjaccard(l, r),
                     getani(l, r)] +
                     [setsketch_jaccard(l, r, k=k, nb=nb, fss=False, executable=executable) for nb in (8, 4, 2, 1)] +
                     [bbminhash_jaccard(l, r, k=k, nb=nb, fss=True, executable=executable) for nb in (8, 4, 2, 1)], np.float32)

def packed(x):
    l, r, k, sz, executable = x
    return getall(l, r, k=j, sz=sz, executable=executable)

def pargetall(tups, k=17, sz=1024, executable="dashing2", cpu=-1):
    import multiprocessing as mp
    if cpu < 0: cpu = mp.cpu_count()
    with mp.Pool(cpu) as p:
        return np.vstack(p.map(packed, tups.ravel().reshape(-1, 2)))


if __name__ == "__main__":
    if(sys.argv[2:]):
        ap = argparse.ArgumentParser()
        ap.add_argument("table", help="Path to a list of selected genomes with ranged Jaccards")
        ap.add_argument("fnames", help="Path to a list of selected genomes")
        ap.add_argument("-k", type=int, default=17, help="Set k for experiment. If exact representations aren't cached for this value of k, it may take a very long time to run")
        ap.add_argument("-s", type=int, default=10, help="Set start sketch size in log2.")
        ap.add_argument("-S", type=int, default=14, help="Set start sketch size in log2.")
        ap.add_argument("-T", type=int, default=2, help="Set step size for sketch size.")
        args = ap.parse_args()
        res = parsedata(args.table, args.fnames)
        bkts, rngs, tups, total_ids, full_genome_ids = res
        selected = [full_genome_ids[i] for i in total_ids]
        print(len(selected), ", frac = %g" % (len(selected) / len(full_genome_ids)), file=sys.stderr)
        s = ""
        print(np.min(tups), np.max(tups))
        print(tups.shape)
        for tup in tups:
            # print(tup.shape)
            s += ",".join("%s-%s" % (x[0], x[1]) for x in tup)
            s += '\n'
        fullmat = np.vstack([np.vstack([getall(l, r, k=args.k, sz=sz) for l, r in map(lambda x: (full_genome_ids[x[0]], full_genome_ids[x[1]]), tups.reshape(-1, 2))]) for sz in map(lambda x: 1<< x, range(args.s, args.S, args.T))])
        settings = [(l, r, sz) for l, r in tups.reshape(-1, 2) for sz in map(lambda x: 1<< x, range(args.s, args.S, args.T))]
        fullmat.astype(np.float32).tofile("fullmat.f32.%s" % (fullmat.shape))
        with open("settings.txt", "w") as f:
            for st in settings:
                print("%s:%d\t%s:d\t%d\n" % (full_genome_ids[st[0]], st[0], full_genome_ids[st[1]], st[1], st[2]), file=f)
    else:
        print("Running tests, not running experiment", file=sys.stderr)
        parse_bf(sys.argv[1] if sys.argv[1:] else "selected_buckets_100.txt")
