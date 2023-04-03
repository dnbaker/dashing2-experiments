#!/usr/bin/env python3
import sys
import itertools
from jir import getall, parsedata
from argparse import ArgumentParser as AP
ap = AP()
ap.add_argument("lhs")
ap.add_argument("rhs")
ap.add_argument("size", type=int)
ap.add_argument("k", type=int)
ap.add_argument("--executable", default="dashing2")
ap = ap.parse_args()
print(f"{ap.lhs}\t{ap.rhs}\t" + "\t".join(map(str, getall(ap.lhs, ap.rhs, k=ap.k, size=ap.size, executable=ap.executable))))
