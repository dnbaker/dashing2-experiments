import numpy as np
import sys
import argparse
import os
import tempfile
from subprocess import PIPE, check_call
from time import time


use_paper_columns = True
include_nibbles = False
include_bmh = False
include_d2minhash = False
include_sourmash = False
include_soumash_abund = False


def check_output(x):
    from subprocess import Popen, PIPE
    t = 0
    while 1:
        process = Popen(x, shell=True, stdout=PIPE, stderr=PIPE, executable="/bin/bash")
        out, err = process.communicate()
        if process.returncode:
            t += 1
            print("Failed call time #%d '%s', returncode = %d, error = %s" % (t, x, process.returncode, err.decode()),
                  file=sys.stderr)
            if t == 10:
                raise RuntimeError("Failed 10 times to run command " + x)
        else:
            return out


def parse_bf(path, ngenomes):
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
    for outer_tup in tups:
        for inner_tup in outer_tup:
            x, y = inner_tup
            assert x < ngenomes
            assert y < ngenomes
    bkts = np.array(bkts)
    ml = max(map(len, tups))
    total_ids = np.array(sorted(set(cfi(cfi(tups)))))
    for i in range(len(tups)):
        if len(tups[i]) < ml:
            diff = ml - len(tups[i])
            tups[i] += [[0, 0]] * diff
        assert len(tups[i]) == ml
        if len(rngs[i]) < 2:
            rngs[i].append(rngs[i][0])
    assert len(set(map(len, rngs))) == 1
    rngs = np.array(rngs)
    # print(f"Total {len(total_ids)} ids parsed", file=sys.stderr)
    return list(map(np.array, (bkts, rngs, tups, total_ids)))


def parsedata(path, fn):
    genome_ids = [x.strip() for x in open(fn)]
    bkts, rngs, tups, ids = parse_bf(path, len(genome_ids))
    return bkts, rngs, tups, ids, genome_ids


def getani(dry, left, r, ex="fastANI"):
    """For genomes left and r, which must exist as files, estimate ANI using
       fastANI returns a numpy array with [ANI, num, denom] as values.
    """
    cmd = f"{ex} --minFraction 0. -q {left} -r {r} -o /dev/stdout"
    if dry:
        return cmd
    out = check_output(cmd).decode().strip()
    try:
        return np.fromiter(out.split("\t")[-3:], dtype=np.float64)[0]
    except Exception as se:
        return 0.


def getmashji(dry, left, r, *, k, size_in_bits):
    if k > 32:
        return -1.  # Mash doesn't support long kmers
    s_arg = size_in_bits//(32 if k <= 16 else 64)
    cmd = f"mash dist -s {s_arg} -k {k} {left} {r}"
    if dry:
        return cmd
    out = check_output(cmd).decode().strip().split("\n")[-1].split()[-1]
    num, denom = map(float, out.split("/"))
    if not denom:
        return 0.
    return num / denom


def getsourmashji(dry, left, r, *, k, size_in_bits):
    s_arg = size_in_bits//64
    with tempfile.TemporaryDirectory() as dr:
        cmd = f"sourmash sketch dna -p k={k},noabund,num={s_arg} {left} {r} -o {dr}/tmp.sig ; sourmash compare {dr}/tmp.sig"
        if dry:
            return cmd
        out = check_output(cmd).decode().strip().split("\n")[-1].split()[-1]
        dr.cleanup()
    return float(out)


def getdashingji(dry, left, r, *, k, l2s=10):
    cyclic_flag = "" if k <= 32 else " --use-cyclic-hash "
    cmd = f"dashing dist {cyclic_flag} -S{l2s} -k{k} {left} {r}"
    if dry:
        return cmd
    return float(check_output(cmd).decode().strip().split("\n")[-2].split("\t")[-1])


def exact_wjaccard(dry, p1, p2, *, k):
    cmd = f"dashing2 sketch --phylip -k{k} --countdict --cache --cmpout /dev/stdout {p1} {p2}"
    if dry:
        return cmd
    return float(check_output(cmd).decode().strip('\n').split('\n')[-2].split("\t")[1])


def exact_jaccard(dry, p1, p2, *, k):
    cmd = f"dashing2 sketch --phylip -k{k} --set --cache --cmpout /dev/stdout {p1} {p2}"
    if dry:
        return cmd
    res = check_output(cmd).decode().strip('\n')
    rf = float(res.split('\n')[-2].split("\t")[1])
    if rf > 1.:
        print("Res", res, rf, file=sys.stderr)
    return rf


def fsarg2str(x):
    if x:
        return " --full-setsketch "
    else:
        return " --oneperm "


def setsketch_jaccard(dry, p1, p2, size_in_bits, *, k, nb=8, fss=False, executable="dashing2"):
    """k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       nb is 8 by default, which is Dashing 2's typical size, though dashing-f and dashing-ld use 4 bytes and 16 bytes
       to store float32 and long double hash registers, respectively.
    """
    # NOTE: nb is expressed in bytes, not bits
    nb_in_bits = nb*8
    cmd = f"{executable} sketch --phylip -k{k} {fsarg2str(fss)} -S{size_in_bits // nb_in_bits} --fastcmp {nb} --cache --cmpout /dev/stdout {p1} {p2}"
    if dry:
        return cmd
    else:
        return float(check_output(cmd).decode().strip('\n').split('\n')[-2].split("\t")[1])


def bbminhash_jaccard(dry, p1, p2, size, *, k, nb=32, fss=False, executable="dashing2"):
    """k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       nb is 8 by default, which is Dashing 2's typical size
    """
    assert (nb & (nb - 1)) == 0 and 4 <= nb <= 64, "nb, number of bits for bbit minhash, should be 4, 8, 16, 32, or 64"
    fcstr = f"--fastcmp {nb / 8}"
    cmd = f"{executable} sketch --phylip -k{k} {fsarg2str(fss)} --bbit-sigs -S{size} {fcstr} --cache --cmpout /dev/stdout {p1} {p2}"
    if dry:
        return cmd
    else:
        return float(check_output(cmd).decode().strip('\n').split('\n')[-2].split("\t")[1])


def bindash_jaccard(dry, p1, p2, size_in_bits, *, k, nb=8, executable="bindash"):
    """Get bindash Jaccard, as well as caching its sketches
    """
    assert size_in_bits % 64 == 0, f"Size must be divisible by 64, found {size_in_bits}"
    bb = nb * 8
    if not os.path.isdir("bdsh_tmp"):
        os.mkdir("bdsh_tmp")

    def mc(p):
        return "bdsh_tmp/" + os.path.basename(p + f".{k}.{bb}.{size_in_bits}")

    cp1, cp2 = map(mc, (p1, p2))
    dry_ret = []
    for cp, p in zip((cp1, cp2), (p1, p2)):
        cmd = f"{executable} sketch --kmerlen={k} --bbits={bb} --minhashtype=2 --sketchsize64={size_in_bits // (64 * bb)} --outfname={cp} {p}"
        if dry:
            dry_ret.append(cmd)
        else:
            if not all(map(os.path.isfile, (cp, cp + ".dat", cp + ".txt"))):
                check_call(cmd, shell=True, stderr=PIPE)
    try:
        cmd = f"{executable} dist {cp1} {cp2}"
        if dry:
            dry_ret.append(cmd)
        else:
            out = check_output(cmd).decode().strip()
    except RuntimeError: # Failed multiple times
        check_call(f"rm {cp1}* {cp2}*", shell=True)
        cmd = f"{executable} sketch --kmerlen={k} --bbits={bb} --minhashtype=2 --sketchsize64={size_in_bits // 64} --outfname={cp} {p}"
        check_call(cmd, shell=True, stderr=PIPE)
        out = check_output(f"{executable} dist {cp1} {cp2}").decode().strip()
    if dry:
        return '\n'.join(dry_ret)
    # print("Output: ", out)
    if not out:
        cmd = f"{executable} dist --ithres=1 {cp1} {cp2}"
        return -1
    toks = out.split('\t')[-1].split("/")
    num, denom = map(float, toks)
    return num / denom


def bagminhash_jaccard(dry, p1, p2, size_in_bits, *, k, nb=8, cssize=-1):
    css = "--countsketch-size %d" % cssize if cssize > 0 else ""
    fss = f" --fastcmp {nb}" if nb in (8, 4, 2, 1, .5) else ""
    cstr = f"dashing2 sketch --phylip --cache --bbit-sigs --cmpout /dev/stdout -S {size_in_bits // (nb * 8)} -k {k} --multiset {css + fss} {p1} {p2}"
    if dry:
        return cstr
    else:
        return float(check_output(cstr).decode().strip('\n').split('\n')[-2].split("\t")[1])


def probminhash_jaccard(dry, p1, p2, size_in_bits, *, k, nb=8, cssize=-1):
    css = "--countsketch-size %d" % cssize if cssize > 0 else ""
    fss = f" --fastcmp {nb}" if nb in (8, 4, 2, 1, .5) else ""
    cstr = f"dashing2 sketch --phylip --cache --bbit-sigs --cmpout /dev/stdout -S {size_in_bits // (nb * 8)} -k {k} --prob {css + fss} {p1} {p2}"
    if dry:
        return cstr
    else:
        return float(check_output(cstr).decode().strip('\n').split('\n')[-2].split("\t")[1])


columns = ['G1', 'G2', 'K', 'sketchsize', 'ANI', 'WJI', 'JI', 'Mash', 'Dash1']
columns = ['G1', 'G2', 'K', 'sketchsize', 'ANI', 'WJI', 'JI', 'Mash', 'Dash1']
if use_paper_columns:
    columns += ['BD1', 'FSS1', 'FSS8', 'SS1', 'SS8', 'PMH1-50000000', 'SM']
else:
    columns += ['BD8', 'BD4', 'BD2', 'BD1']
    if include_nibbles:
        columns += ['BDN']
    columns += ['SS8', 'SS4', 'SS2', 'SS1']
    if include_nibbles:
        columns += ['SSN']
    columns += ['FSS8', 'FSS4', 'FSS2', 'FSS1']
    if include_nibbles:
        columns += ['FSSN']
    if include_d2minhash:
        columns += ['MH8', 'MH4', 'MH2', 'MH1']
        if include_nibbles:
            columns += ['MHN']
        columns += ['FMH8', 'FMH4', 'FMH2', 'FMH1']
        if include_nibbles:
            columns += ['FMHN']

    if include_sourmash:
        columns += ['SM']
    if include_soumash_abund:
        columns += ['SMA']

    # In the paper we only look at examples where the CountMinSketch has 5M 64-bit counts
    PMNBs = [8, 4, 2, 1]
    BMNBs = [8, 4, 2, 1]
    if include_nibbles:
        PMNBs.append(.5)
        BMNBs.append(.5)
    CSSZ = [-1, 50000000, 500000]
    PMHSettings = [(b, cs) for b in PMNBs for cs in CSSZ]
    BMHSettings = [(b, cs) for b in BMNBs for cs in CSSZ]
    for (b, cs) in PMHSettings:
        columns.append("PMH%s%s" % (b if b >= 1 else "N", "-%d" % cs if cs > 0 else "Exact"))
    if include_bmh:
        for (b, cs) in BMHSettings:
            columns.append("BMH%s%s" % (b if b >= 1 else "N", "-%d" % cs if cs > 0 else "Exact"))

header = "#" + '\t'.join(columns)
ncols = len(columns)

bbnbs = [8, 4, 2, 1]
if include_nibbles:
    bbnbs.append(.5)


def getall(dry, l, r, k=17, size_in_bits=1024, executable="dashing2", faex="fastANI"):
    """
        for values of k, size, and executable, return
        all similarity comparisons using Mash, Dashing, Dashing2, and fastANI
    """
    results = []
    if 'ANI' in columns:
        results += [getani(dry, l, r, ex=faex)]
    if 'WJI' in columns:
        results += [exact_wjaccard(dry, l, r, k=k)]
    if 'JI' in columns:
        results += [exact_jaccard(dry, l, r, k=k)]
    if 'Mash' in columns:
        results += [getmashji(dry, l, r, k=k, size_in_bits=size_in_bits)]
    if 'Dash1' in columns:
        results += [getdashingji(dry, l, r, k=k, l2s=int(np.log2(size_in_bits//8)))]
    if 'SM' in columns:
        results += [getsourmashji(dry, l, r, k=k, size_in_bits=size_in_bits)]
    if 'BD1' in columns:
        results += [bindash_jaccard(dry, l, r, k=k, size_in_bits=size_in_bits, nb=nb) for nb in bbnbs]
    if 'SS1' in columns:
        results += [setsketch_jaccard(dry, l, r, size_in_bits=size_in_bits, k=k, nb=nb, fss=False, executable=executable) for nb in bbnbs]
    if 'FSS1' in columns:
        results += [setsketch_jaccard(dry, l, r, size_in_bits=size_in_bits, k=k, nb=nb, fss=True, executable=executable) for nb in bbnbs]
    if 'MH1' in columns:
        results += [bbminhash_jaccard(dry, l, r, size=size, k=k, nb=int(nb * 8), fss=False, executable=executable) for nb in bbnbs]
    if 'FMH1' in columns:
        results += [bbminhash_jaccard(dry, l, r, size=size, k=k, nb=int(nb * 8), fss=True, executable=executable) for nb in bbnbs]
    if 'PMH1Exact' in columns:
        results += [probminhash_jaccard(dry, l, r, size_in_bits=size_in_bits, k=k, nb=b, cssize=csz) for b, csz in PMHSettings]
    if 'BMH1Exact' in columns:
        results += [bagminhash_jaccard(dry, l, r, size_in_bits=size_in_bits, k=k, nb=b, cssize=csz) for b, csz in BMHSettings]
    ret = np.array(results, np.float32)
    print(ret)
    if not dry and ret[2] > 1.:
        print(f"Distance {l} {r}, k = {k}, size = {size}... For some reason, Jaccard was > 1...What's going on?",
              ret[2], file=sys.stderr)
    return ret


def packed(x):
    return packed_dry(False, x)


def packed_dry(dry, x):
    try:
        l, r, k, size_in_bits, executable = x
    except Exception as _:
        print(len(x), file=sys.stderr)
        raise
    ret = getall(dry, l, r, k=k, size_in_bits=size_in_bits, executable=executable)
    return ret


if __name__ == "__main__":
    if sys.argv[2:] or "-h" in sys.argv or "--help" in sys.argv:
        ap = argparse.ArgumentParser()
        ap.add_argument("table", help="Path to a list of selected genomes with ranged Jaccards")
        ap.add_argument("fnames", help="Path to a list of genomes identifiers")
        ap.add_argument("-k", type=int, default=17, help="Set k for experiment. If exact representations aren't cached for this value of k, it may take a very long time to run")
        ap.add_argument("-s", type=int, default=15, help="Set start sketch size in log2 bits.")
        ap.add_argument("-S", type=int, default=19, help="Set end sketch size in log2 bits.")
        ap.add_argument("-T", type=int, default=1, help="Set step size for sketch size.")
        ap.add_argument("--cpu", type=int, default=-1)
        ap.add_argument("--executable", '-E', default="dashing2")
        ap.add_argument("--name", default="noname")
        ap.add_argument("--dryrun", action='store_true')
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
        if not args.dryrun:
            rawmat = []
            with mp.Pool(cpu) as p:
                for i, st in enumerate(subtups):
                    startt = time()
                    print("Started subgroup %d/%d" % (i, len(subtups)), file=sys.stderr, flush=True)
                    result_list = list(p.map(packed, st))
                    for (left, r, k, size, _), mr in zip(st, result_list):
                        print(f"{left}\t{r}\t{k}\t{size}\t" + "\t".join(map(str, mr)), file=ofp, flush=True)
                    rawmat.append(np.stack(result_list))
            if ofp is not sys.stdout:
                ofp.close()
            np.vstack(rawmat).tofile(f"fullmat.{args.name}.f32.{hv}.{k}.{fs}")
        else:
            for st in subtups:
                print(list(map(lambda x: packed_dry(True, x), st)))
    else:
        print("Running tests, not running experiment", file=sys.stderr)
        parse_bf(sys.argv[1] if sys.argv[1:] else "selected_buckets_100.txt")
