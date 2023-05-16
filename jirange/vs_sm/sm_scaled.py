#!/usr/bin/env python3

"""
Script for running
"""

import os
import sys
import json
import time
import numpy as np
import argparse


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
        print('Skipping creation of "%s" as it already exists' % sig_file)
    with open(sig_file, 'r') as fh:
        return len(json.loads(fh.read())[0]['signatures'][0]['mins'])


def sourmash_compare(in_prefix, k=21, force=False):
    np_file = os.path.join(in_prefix, 'sm_ani_containment.np')
    in1 = os.path.join(in_prefix, '1.fna.gz')
    in2 = os.path.join(in_prefix, '2.fna.gz')
    sig_suffix = '.' + '.'.join(['k%d' % k, 'sig'])
    assert os.path.exists(in1 + sig_suffix)
    assert os.path.exists(in2 + sig_suffix)
    if not os.path.exists(np_file) or force:
        cmd = ['sourmash', 'compare', '--ani', '--containment', in1 + sig_suffix, in2 + sig_suffix, '-o', np_file]
        ret = os.system(' '.join(cmd))
        if ret != 0:
            raise RuntimeError('Command "%s" returned non-zero %d' % (' '.join(cmd), ret))
        if not os.path.exists(np_file):
            raise RuntimeError('Numpy file "%s" does not exist' % np_file)
    return np.load(np_file)[0, 1]


def check_output(x, retries=3):
    from subprocess import Popen, PIPE
    t = 0
    while 1:
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


def fastani(in_prefix, ex="fastANI", force=False):
    """For genomes left and r, which must exist as files, estimate ANI using
       fastANI returns a numpy array with [ANI, num, denom] as values.
    """
    in1 = os.path.join(in_prefix, '1.fna.gz')
    in2 = os.path.join(in_prefix, '2.fna.gz')
    fani_file = os.path.join(in_prefix, '.fastani')
    if not os.path.exists(fani_file) or force:
        cmd = f"{ex} --minFraction 0. -q {in1} -r {in2} -o {fani_file}"
        ret = os.system(' '.join(cmd))
        if ret != 0:
            raise RuntimeError('Command "%s" returned non-zero %d' % (' '.join(cmd), ret))
        if not os.path.exists(fani_file):
            raise RuntimeError('FastANI output file "%s" does not exist' % fani_file)
    with open(fani_file, 'rt') as fh:
        ln = fh.readline().rstrip()
        return ln.split('\t')


def dashing2_sketch(inf, k=21, force=False):
    d2_file = '.'.join([inf, 'k%d' % k, 'op'])
    if not os.path.exists(d2_file) or force:
        sketch - -binary - output - S
        20 - k
        31
        cmd = f"dashing2 sketch --phylip -k{k} --set --cache --cmpout /dev/stdout {p1} {p2}"


def dashing2_compare():
    in1 = os.path.join(in_prefix, '1.fna.gz')
    in2 = os.path.join(in_prefix, '2.fna.gz')
    fani_file = os.path.join(in_prefix, '.fastani')
    cmd = f"dashing2 sketch --phylip -k{k} --set --cache --cmpout /dev/stdout {p1} {p2}"


def go():
    if not os.path.exists('pairs'):
        raise RuntimeError('No "pairs" subdirectory')
    for i in range(1, 1011):
        istr = '%04d' % i
        subdir = os.path.join('pairs', istr)
        in1 = os.path.join(subdir, '1.fna.gz')
        in2 = os.path.join(subdir, '2.fna.gz')
        if not os.path.exists(in1) or not os.path.exists(in2):
            raise RuntimeError('failed to find an expected file "%s" or "%s"' % (in1, in2))
        sk1_sz = sourmash_sketch(in1)
        sk2_sz = sourmash_sketch(in2)
        jaccard = sourmash_compare(subdir)
        fastani = fastani(in1)
        print("sk1sz=%d, sk2sz=%d, J=%f" % (sk1_sz, sk2_sz, jaccard))


if __name__ == '__main__':
    go()
