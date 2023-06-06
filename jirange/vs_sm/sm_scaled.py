#!/usr/bin/env python3

"""
Script for running
"""

import os
import sys
import json
import time
import math
import glob
import numpy as np
import multiprocessing


# df$ani_est <- ifelse(df$jest > 0, 1 + 1/df$k * log(2*df$jest/(1+df$jest)), 0)
# df$ani_est <- pmax(df$ani_est, 0.0)
def jaccard_to_ani(j, k):
    if j <= 0:
        return 0
    return max(0, 1 + (1/k) * math.log((2 * j) / (1 + j)))


def sourmash_sketch(inf, k=21, force=False):
    sig_file = '.'.join([inf, 'k%d' % k, 'sig'])
    if not os.path.exists(sig_file) or force:
        cmd = ['sourmash', 'sketch', 'dna', '-p', 'k=%d' % k, inf, '-o', sig_file]
        ret = os.system(' '.join(cmd))
        if ret != 0:
            raise RuntimeError('Command "%s" returned non-zero %d' % (' '.join(cmd), ret))
        if not os.path.exists(sig_file):
            raise RuntimeError('Signature file "%s" does not exist' % sig_file)
    else:
        print('Skipping creation of "%s" as it already exists' % sig_file, file=sys.stderr)
    with open(sig_file, 'r') as fh:
        return len(json.loads(fh.read())[0]['signatures'][0]['mins'])


def sourmash_compare(in_prefix, k=21, force=False):
    np_file = os.path.join(in_prefix, f'sm_ani_containment.k{k}.np')
    in1 = os.path.join(in_prefix, '1.fna.gz')
    in2 = os.path.join(in_prefix, '2.fna.gz')
    sig_suffix = '.' + '.'.join(['k%d' % k, 'sig'])
    assert os.path.exists(in1 + sig_suffix)
    assert os.path.exists(in2 + sig_suffix)
    if not os.path.exists(np_file) or force:
        cmd = ['sourmash', 'compare', '--ani', '--containment', in1 + sig_suffix, in2 + sig_suffix, '-o', np_file]
        ret = os.system(' '.join(cmd) + ' 2>/dev/null')
        if ret != 0:
            raise RuntimeError('Command "%s" returned non-zero %d' % (' '.join(cmd), ret))
        if not os.path.exists(np_file):
            raise RuntimeError('Numpy file "%s" does not exist' % np_file)
    return np.load(np_file)[0, 1]


def check_output(x, retries=3, verbose=True):
    from subprocess import Popen, PIPE
    t = 0
    while 1:
        if verbose:
            print(x)
        process = Popen(x, shell=True, stdout=PIPE, stderr=PIPE, executable="/bin/bash")
        out, err = process.communicate()
        if process.returncode:
            t += 1
            print("Failed call time #%d '%s', returncode = %d, error = %s" %
                  (t, x, process.returncode, err.decode()),
                  file=sys.stderr)
            time.sleep(1)
            if t == retries:
                raise RuntimeError("Failed 10 times to run command " + x)
        else:
            return out


def fastani(in_prefix, ex="fastANI", force=False, verbose=True):
    """For genomes left and r, which must exist as files, estimate ANI using
       fastANI returns a numpy array with [ANI, num, denom] as values.
    """
    in1 = os.path.join(in_prefix, '1.fna.gz')
    in2 = os.path.join(in_prefix, '2.fna.gz')
    fani_file = os.path.join(in_prefix, 'fastani.tsv')
    if not os.path.exists(fani_file) or force:
        cmd = f"{ex} --minFraction 0. -q {in1} -r {in2} -o {fani_file} 2>/dev/null"
        if verbose:
            print(cmd)
        ret = os.system(cmd)
        if ret != 0:
            raise RuntimeError('Command "%s" returned non-zero %d' % (' '.join(cmd), ret))
        if not os.path.exists(fani_file):
            raise RuntimeError('FastANI output file "%s" does not exist' % fani_file)
    with open(fani_file, 'rt') as fh:
        ln = fh.readline().rstrip()
        return float(ln.split('\t')[-3])


def fsarg2str(x):
    if x:
        return " --full-setsketch "
    else:
        return " --oneperm "


def dashing2_sketch(inf, size_in_bytes, k=21, full=False, cssize=-1, regsize=1):
    css = "--prob --countsketch-size %d" % cssize if cssize > 0 else ""
    cmd = f"dashing2 sketch {css} -k{k} {fsarg2str(full)} -S{size_in_bytes} " + \
          f"--cache --fastcmp {regsize} {inf} 2>/dev/null"
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f'Nonzero exit {ret} from "{cmd}"')


def dashing2_compare(in_prefix, size_in_bytes, k=21, full=False, regsize=1, cssize=-1):
    in1 = os.path.join(in_prefix, '1.fna.gz')
    in2 = os.path.join(in_prefix, '2.fna.gz')

    flat_sketches = glob.glob('in_prefix/*.opss')
    full_sketches = glob.glob('in_prefix/*.ss')
    weighted_sketches = glob.glob('in_prefix/*.pmh')

    css = "--prob --countsketch-size %d" % cssize if cssize > 0 else ""
    cmd = f"dashing2 cmp --cache -k{k} {fsarg2str(full)} -S{size_in_bytes} " + \
          f"--fastcmp {regsize} {css} --cmpout /dev/stdout {in1} {in2} 2>/dev/null"
    out = check_output(cmd, verbose=False)

    assert len(glob.glob('in_prefix/*.opss')) == len(flat_sketches)
    assert len(glob.glob('in_prefix/*.ss')) == len(full_sketches)
    assert len(glob.glob('in_prefix/*.pmh')) == len(weighted_sketches)

    lines = out.split(b'\n')
    assert len(lines) >= 3
    assert len(lines[-3]) > 0
    j = float(lines[-3].split()[-1])
    return j


def sm_sketch_task(tup):
    _k, _in1, _in2, _subdir = tup
    _sk1_sz = sourmash_sketch(_in1, k=_k)
    _sk2_sz = sourmash_sketch(_in2, k=_k)
    with open(os.path.join(_subdir, f'sm_sk_sz.k{_k}.csv'), 'wt') as _fh:
        _fh.write(f'{_sk1_sz},{_sk2_sz}\n')


def sm_compare_task(tup):
    _k, _, _, _subdir = tup
    return sourmash_compare(_subdir, k=_k)


def d2_sketch_task(tup):
    _k, _in1, _in2, _subdir = tup
    with open(os.path.join(_subdir, f'sm_sk_sz.k{_k}.csv'), 'rt') as _fh:
        _sk1_sz, _sk2_sz = _fh.read().rstrip().split(',')
        _sk1_sz, _sk2_sz = int(_sk1_sz), int(_sk2_sz)
    sk_sz = min(_sk1_sz, _sk2_sz)

    flat = ['rc_canon', f'sketchsize{sk_sz}', f'k{_k}', 'SetSpace', 'DNA', 'opss']
    _fn = '.'.join([_in1] + flat)
    if os.path.exists(_fn):
        print(f'Skipping creating of "{_fn}" as it already exists...', file=sys.stderr)
    else:
        dashing2_sketch(_in1, k=_k, size_in_bytes=sk_sz)

    _fn = '.'.join([_in2] + flat)
    if os.path.exists(_fn):
        print(f'Skipping creating of "{_fn}" as it already exists...', file=sys.stderr)
    else:
        dashing2_sketch(_in2, k=_k, size_in_bytes=sk_sz)

    full = ['rc_canon', f'sketchsize{sk_sz}', f'k{_k}', 'SetSpace', 'DNA', 'ss']
    _fn = '.'.join([_in1] + full)
    if os.path.exists(_fn):
        print(f'Skipping creating of "{_fn}" as it already exists...', file=sys.stderr)
    else:
        dashing2_sketch(_in1, k=_k, full=True, size_in_bytes=sk_sz)

    _fn = '.'.join([_in2] + full)
    if os.path.exists(_fn):
        print(f'Skipping creating of "{_fn}" as it already exists...', file=sys.stderr)
    else:
        dashing2_sketch(_in2, k=_k, full=True, size_in_bytes=sk_sz)

    weighted = ['rc_canon', f'sketchsize{sk_sz}', f'k{_k}', 'CountMinCounting50000000', 'ProbsetSpace', 'DNA', 'pmh']
    _fn = '.'.join([_in1] + weighted)
    if os.path.exists(_fn):
        print(f'Skipping creating of "{_fn}" as it already exists...', file=sys.stderr)
    else:
        dashing2_sketch(_in1, k=_k, cssize=50000000, size_in_bytes=sk_sz)

    _fn = '.'.join([_in2] + weighted)
    if os.path.exists(_fn):
        print(f'Skipping creating of "{_fn}" as it already exists...', file=sys.stderr)
    else:
        dashing2_sketch(_in2, k=_k, cssize=50000000, size_in_bytes=sk_sz)


def d2_compare_task(tup):
    _k, _, _, _subdir = tup
    with open(os.path.join(_subdir, f'sm_sk_sz.k{_k}.csv'), 'rt') as _fh:
        _sk1_sz, _sk2_sz = _fh.read().rstrip().split(',')
        _sk1_sz, _sk2_sz = int(_sk1_sz), int(_sk2_sz)
    sk_sz = min(_sk1_sz, _sk2_sz)
    d2_file = os.path.join(_subdir, f'd2.k{_k}.csv')
    if not os.path.exists(d2_file):
        j1 = dashing2_compare(_subdir, sk_sz, k=_k, full=False, regsize=1, cssize=-1)
        j2 = dashing2_compare(_subdir, sk_sz, k=_k, full=True, regsize=1, cssize=-1)
        j3 = dashing2_compare(_subdir, sk_sz, k=_k, full=True, regsize=1, cssize=50000000)
        a1, a2, a3 = jaccard_to_ani(j1, k=_k), jaccard_to_ani(j2, k=_k), jaccard_to_ani(j3, k=_k)
        with open(d2_file, 'wt') as fh:
            fh.write(','.join(map(str, [a1, a2, a3])))
        return a1, a2, a3
    else:
        with open(d2_file, 'rt') as fh:
            results = list(map(float, fh.read().rstrip().split(',')))
        assert 3 == len(results)
        return results[0], results[1], results[2]


def fastani_task(tup):
    _, _, _, _subdir = tup
    return fastani(_subdir) / 100


def go():
    if not os.path.exists('pairs'):
        raise RuntimeError('No "pairs" subdirectory')

    tasks = []
    for k in [21, 31]:
        for i in range(1, 1011):
            istr = '%04d' % i
            subdir = os.path.join('pairs', istr)
            in1 = os.path.join(subdir, '1.fna.gz')
            in2 = os.path.join(subdir, '2.fna.gz')
            if not os.path.exists(in1) or not os.path.exists(in2):
                raise RuntimeError('failed to find an expected file "%s" or "%s"' % (in1, in2))
            tasks.append((k, in1, in2, subdir))

    print('Starting sourmash SKETCH tasks...', file=sys.stderr)
    with multiprocessing.Pool() as pool:
        pool.map(sm_sketch_task, tasks)
    print('DONE sourmash sketch tasks', file=sys.stderr)

    print('Starting d2 SKETCH tasks...', file=sys.stderr)
    with multiprocessing.Pool() as pool:
        pool.map(d2_sketch_task, tasks)
    print('DONE d2 sketch tasks', file=sys.stderr)

    if not os.path.exists('fastani_results.txt'):
        print('Starting fastANI tasks...', file=sys.stderr)
        with multiprocessing.Pool() as pool:
            fastani_results = pool.map(fastani_task, tasks)
        print(f'DONE {len(fastani_results)} fastANI tasks', file=sys.stderr)
        with open('fastani_results.txt', 'wt') as fh:
            print(fastani_results, file=fh)
    else:
        print('Skipping fastANI COMPARE tasks since fastani_results.txt exists', file=sys.stderr)

    if not os.path.exists('sm_results.txt'):
        print('Starting sourmash COMPARE tasks...', file=sys.stderr)
        with multiprocessing.Pool() as pool:
            sm_results = pool.map(sm_compare_task, tasks)
        print(f'DONE {len(sm_results)} sourmash compare tasks', file=sys.stderr)
        with open('sm_results.txt', 'wt') as fh:
            print(sm_results, file=fh)
    else:
        print('Skipping sourmash COMPARE tasks since sm_results.txt exists', file=sys.stderr)

    if not os.path.exists('d2_results.txt'):
        print('Starting d2 COMPARE tasks...', file=sys.stderr)
        with multiprocessing.Pool() as pool:
            d2_results = pool.map(d2_compare_task, tasks)
        print(f'DONE {len(d2_results)} d2 compare tasks', file=sys.stderr)
        with open('d2_results.txt', 'wt') as fh:
            print(d2_results, file=fh)
    else:
        print('Skipping Dashing 2 COMPARE tasks since d2_results.txt exists', file=sys.stderr)

    print(','.join(['k', 'sm_ani', 'd2_ani', 'd2f_ani', 'd2w_ani', 'fastani_ani']))
    for tup in zip(tasks, fastani_results, sm_results, d2_results):
        task_tup, fastani_result, sm_result, d2_result_tup = tup
        d2_result, d2f_result, d2w_result = d2_result_tup
        k, _, _, _ = task_tup
        print(','.join(map(str, [k, sm_result, d2_result, d2f_result, d2w_result, fastani_result])))


if __name__ == '__main__':
    go()
