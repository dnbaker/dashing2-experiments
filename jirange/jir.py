import numpy as np
import scipy.sparse as sp
import sys
import subprocess


# 1. Parse the ids
# 2. For each relevant id, perform full k-mer count dictionary generation
# 3. Then, we can use this to calculate shared k-mers, weighted k-mer comparisons.
# 4. Then, we calculate all-pairs over the relevant genomes for all of these.
#   1. This includes exact distances, the range of setsketch parameters, HLL, Mash, Dashing1, bindash, ARI
# 5. Then, parse those results out and prepare the table or results

# This will give us ANI, jaccard, etc. and will do super awesome powerful stuff

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
    print(len(selected_genomes))
    return (bkts, rngs, tups, ids, genome_ids)


def getani(l, r): 
    '''For genomes l and r, which must exist as files, estimate ANI using fastANI
       returns a numpy array with [ANI, num, denom] as values.
    '''
    from subprocess import check_output
    # return [ANI, num, denom]
    cmd = f"fastANI --minFraction 0. -q {l} -r {r} -o /dev/stdout"
    return np.fromiter(check_output(cmd, shell=True).decode().strip().split("\t")[-3:], dtype=np.float64)


def getmashji(l, r, *, k):
    return float(check_output(f"jaccard_mash dist -k {k} -j -t {l} {r}", shell=True).decode().strip().split("\n")[1].split()[-1])


# To add: jaccard_mash (the mash fork which has the -j CLI argument)


def exact_wjaccard(p1, p2, k=17):
    return float(check_output("dashing2 sketch -k{k} --countdict --exact-kmer-dist --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])

def exact_jaccard(p1, p2, k=17):
    return float(check_output("dashing2 sketch -k{k} --set --exact-kmer-dist --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def setsketch_jaccard(p1, p2, size, k=17, nb=8, fss=False):
    '''k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       nb is 8 by default, which is Dashing 2's typical size, though dashing-f and dashing-ld use 4 bytes and 16 bytes
       to store float32 and long double hash registers, respectively.
    '''
    return float(check_output("dashing2 sketch -k{k} {'--full-setsketch ' if fss else ''} -S{size} --fastcmp {nb} --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])

def bbminhash_jaccard(p1, p2, size, k=17, bbnb=32, fss=False):
    '''k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       bbnb is 8 by default, which is Dashing 2's typical size
    '''
    nb = bbnb // 8
    assert (bbnb & (bbnb - 1)) == 0 and 4 <= bbnb <= 64, "bbnb, number of bits for bbit minhash, should be 4, 8, 16, 32, or 64"
    return float(check_output("dashing2 sketch -k{k} {'--full-setsketch ' if fss else ''} --bbit-sigs -S{size} --fastcmp {nb} --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


if __name__ == "__main__":
    parse_bf(sys.argv[1] if sys.argv[1:] else "selected_buckets_100.txt")
    if(sys.argv[2:]):
        parsedata(sys.argv[1], sys.argv[2])

