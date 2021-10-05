import sys
import os
import numpy as np
import subprocess
from timeit import default_timer as time
import random

def getrstr(x):
    return "".join(map(str, (random.choice(range(10)) for x in range(20))))


def bch_sketch_mash(pathf, k, threads, destp, size):
    starttm = time()
    errfile = getrstr()
    try:
        subprocess.check_call(f"mash sketch -s {size} -k {k} -o {destp} -l {pathf} -p {threads} 2> {errfile}", shell=True)
    except subprocess.CalledProcessError as e:
        print(f"error: {e}. Errfile: {errfile}", file=sys.stderr)
        raise
    stoptm = time()
    destp += ".msh"
    assert os.path.isfile(destp)
    return destp, (stoptm - starttm)


def bch_sketch_dashing(pathf, k, threads, size):
    starttd1 = time()
    errfile = getrstr()
    try:
        subprocess.check_call(f"dashing sketch -S {int(np.ceil(np.log2(size)))} -k {k} -F {pathf} -p {threads}, 2>errfile.{errfile}.txt", shell=True)
    except subprocess.CalledProcessError as e:
        print(f"error: {e}. Errfile: {errfile}", file=sys.stderr)
        raise
    stoptd1 = time()
    return pathf, (stoptd1 - starttd1)


def bch_sketch_dashing2(pathf, k, threads, size, oneperm=False):
    starttd2 = time()
    pstr = " --oph " if oneperm else ""
    errfile = getrstr()
    try:
        subprocess.check_call(f"dashing2 sketch -k {k} -S {size} -F {pathf} -p {threads} {pstr} 2>errfile.{errfile}.txt ", shell=True)
    except subprocess.CalledProcessError as e:
        print(f"error: {e}. Errfile: {errfile}", file=sys.stderr)
        raise
    return pathf, (time() - starttd2)


def bch_dist_dashing2(pathf, k, threads, size, oneperm=False, distdest=None, binary=False, regsize=-1):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    if regsize <= 0. or regsize not in [0.5, 1, 2, 4, 8]:
        raise ValueError("regsize must be > 0 and in [.5, 1, 2, 4, 8]")
    startt = time()
    pstr = " --oph " if oneperm else ""
    pstr += " --binary-output " if binary else ""
    errfile = getrstr()
    try:
        subprocess.check_call(f"dashing2 sketch -k {k} --fastcmp {regsize} --cache --cmpout {distdest} -S {size} -F {pathf} -p {threads} {pstr} 2>errfile.{errfile}.txt", shell=True)
    except subprocess.CalledProcessError as e:
        print(f"error: {e}. Errfile: {errfile}", file=sys.stderr)
        raise
    return pathf, (time() - startt)


def bch_dist_dashing(pathf, k, threads, size, distdest=None, binary=False):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    startt = time()
    pstr = " --emit-binary" if binary else ""
    subprocess.check_call(f"dashing dist --cache-sketches -k {k} -O{distdest} -o {distdest + '.sizes'} -S {int(np.ceil(np.log2(size)))} -F {pathf} -p {threads} {pstr} 2>errfile.{getrstr()}.txt", shell=True)
    return pathf, (time() - startt)


def bch_dist_bindash(pathf, threads, distdest=None):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    startt = time()
    subprocess.check_call(f"bindash dist --outfname={distdest} --nthreads={threads} {pathf} 2>errfile.{getrstr()}.txt", shell=True)
    return pathf, (time() - startt)


def bch_dist_mash(pathf, threads, distdest=None):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_dashing2")
    startt = time()
    subprocess.check_call(f"mash triangle -p {threads} {pathf} > {distdest}", shell=True)
    return pathf, (time() - startt)


def bch_sketch_bindash(pathf, k, threads, *, destp, bbits, size):
    startt = time()
    subprocess.check_call(f"bindash sketch --kmerlen={k} --sketchsize64={size//64} --listfname={pathf} --outfname={destp} --nthreads={threads} 2>errfile.{getrstr()}.txt", shell=True)
    return destp, (time() - startt)


def repeat_x(func, ntimes, *args, **kwargs):
    timevec = np.zeros(ntimes)
    results = [func(*args, **kwargs) for i in range(ntimes)]
    timevec = [x[-1] for x in results]
    return (results[0][0], np.median([x[-1] for x in results]))


def bch_sketch_pmh(pathf, k, threads, size, *, cssize=None):
    starttpmh = time()
    csstr = f"--countsketch-size {cssize}" if cssize is not None else ""
    subprocess.check_call(f"dashing2 sketch --probminhash -k {k} -S {size} -F {pathf} -p {threads} {csstr} 2>errfile.{getrstr()}.txt ", shell=True)
    return pathf, (time() - starttpmh)


def bch_dist_pmh(pathf, k, threads, size, *, cssize=None, distdest=None, binary=False, regsize=-1):
    if distdest is None:
        raise RuntimeError("Must provide distdest to bch_dist_pmh")
    if regsize <= 0. or regsize not in [0.5, 1, 2, 4, 8]:
        raise ValueError("regsize must be > 0 and in [.5, 1, 2, 4, 8]")
    startt = time()
    pstr = " --binary-output " if binary else ""
    pstr += f" --countsketch-size {cssize}" if cssize is not None else ""
    subprocess.check_call(f"dashing2 sketch -k {k} --bbit-sigs --fastcmp {regsize} --cache --cmpout {distdest} -S {size} -F {pathf} -p {threads} {pstr} 2>errfile.{getrstr()}.txt", shell=True)
    return pathf, (time() - startt)
