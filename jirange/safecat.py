#!/usr/bin/env python3
import sys
import os


def extractfromrn(x):
    x = x.strip()
    if not x:
        return None
    toks = x.strip().split('.')
    # assert toks[-1] == "row"
    sketchsize, k = toks[-4:-2]
    with open(x) as ifp:
        try:
            l = next(ifp).strip().split()
        except StopIteration:
            return None
    from itertools import chain
    return "\t".join(chain(l[:2], (k, sketchsize), l[2:]))


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("inputfile", help="Document containing the list of row-files")
    ap.add_argument("outputfile", help="Path to which to write the final output table")
    ap.add_argument("--linenum", default=0, type=int)
    header = "#G1\tG2\tK\tsketchsize\tANI\tWJI\tJI\tMash\tDash1\tBD8\tBD4\tBD2\tBD1\tBDN\tSS8\tSS4\tSS2\tSS1\tSSN\tFSS8\tFSS4\tFSS2\tFSS1\tFSSN\tMH8\tMH4\tMH2\tMH1\tMHN\tFMH8\tFMH4\tFMH2\tFMH1\tFMHN\tPMH8Exact\tPMH8-50000000\tPMH8-500000\tPMH4Exact\tPMH4-50000000\tPMH4-500000\tPMH2Exact\tPMH2-50000000\tPMH2-500000\tPMH1Exact\tPMH1-50000000\tPMH1-500000\tPMHNExact\tPMHN-50000000\tPMHN-500000\tBMH8Exact\tBMH8-50000000\tBMH8-500000\tBMH4Exact\tBMH4-50000000\tBMH4-500000\tBMH2Exact\tBMH2-50000000\tBMH2-500000\tBMH1Exact\tBMH1-50000000\tBMH1-500000\tBMHNExact\tBMHN-50000000\tBMHN-500000"
    ap = ap.parse_args()
    outmode = "w" if ap.linenum == 0 else "a"
    with open(ap.outputfile, outmode) as ofp, open(ap.inputfile) as ifp:
        if ap.linenum == 0:
            print(header, file=ofp)
        for i in range(ap.linenum):
            if (i & 1023) == 0:
                print("Skipped %d/%d lines ahead" % (i, ap.linenum), file=sys.stderr)
            _ = next(ifp) # Skip that many lines
        def print_if(item):
            item = extractfromrn(item)
            if item:
                print(item, file=ofp)
        set(map(print_if, ifp))
