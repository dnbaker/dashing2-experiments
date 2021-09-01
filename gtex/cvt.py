import numpy as np
import scipy.sparse as sp
import argparse as ap
import os
import sys

# Script for converting back

aa = ap.ArgumentParser()
aa.add_argument("prefix", nargs='?', default="GTex", help="Provide input prefix; defaults to GTex. Determines both input and output.")
aa = aa.parse_args()


def getmem():
    from psutil import Process
    return Process(os.getpid()).memory_info().rss

print(f"getmem before checking ctsv ax value: {getmem()}", file=sys.stderr)
ctsv = np.memmap(aa.prefix + ".cts.u32", np.uint32)
is32 = True
if np.max(ctsv) < 65536:
    ctsv = ctsv.astype(np.uint16)
    is32 = False
print(f"getmem before making matrix: {getmem()}", file=sys.stderr)
mat = sp.csr_matrix(ctsv, np.memmap(aa.prefix + ".ids.u16", np.uint16), np.memmap(aa.prefix + ".indptr.u64", np.uint64)).T.tocsr()
mat.sort_indices()
print(f"getmem before final conversions: {getmem()}", file=sys.stderr)
mat.data = mat.data.astype(np.uint32 if is32 else np.uint16)
mat.data.tofile(aa.prefix + ".T.cts.u" + "32" if is32 else "16")
maxind = np.max(mat.indices)
if maxind <= 0xFF:
    mat.indices.astype(np.uint8).tofile(aa.prefix + ".T.ids.u8")
elif maxind <= 0xFFFF:
    mat.indices.astype(np.uint16).tofile(aa.prefix + ".T.ids.u8")
elif maxind <= 0xFFFFFFFF:
    mat.indices.astype(np.uint32).tofile(aa.prefix + ".T.ids.u32")
else:
    mat.indices.astype(np.uint64).tofile(aa.prefix + ".T.ids.u64")
mat.indptr.astype(np.uint64).tofile(aa.prefix + ".T.indptr.u64")
np.array(mat.shape, np.uint64).tofile(aa.prefix + ".T.shape.u64")
print(f"getmem at end: {getmem()}", file=sys.stderr)
