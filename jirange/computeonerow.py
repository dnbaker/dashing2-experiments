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
ap = ap.parse_args()
lhs, rhs = ap.lhs, ap.rhs
size = ap.size
k = ap.k
print(f"{lhs}\t{rhs}\t" + "\t".join(map(str, getall(lhs, rhs, k=k, size=size))))
