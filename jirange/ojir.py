import numpy as np
import scipy.sparse as sp
import sys


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


if __name__ == "__main__":
    parse_bf(sys.argv[1])
    if(sys.argv[2:]):
        parsedata(sys.argv[1], sys.argv[2])
