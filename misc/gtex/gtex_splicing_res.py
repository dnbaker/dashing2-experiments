import scipy.sparse as sp
import minicore as mc
import numpy as np
import lzma
import argparse as ap
ap = ap.ArgumentParser()
aa = ap.add_argument
aa("prior", type=float, default=ap.prior)
ap = ap.parse_args()

mat = sp.csr_matrix((np.fromfile("jids.data.u16.npy", np.uint16), np.fromfile("jids.indices.u32.npy"), np.uint32), np.fromfile("jids.indptr.u64.npy", np.uint64))
print("Loaded matrix")
mat.indices = mat.indices.astype(np.uint32)
mat.indptr = mat.indptr.astype(np.uint64)
mc.set_num_threads(32)
csm = mc.CSparseMatrix(mat)
sub = ["32", "64"][pdists.dtype.char == 'd']
pdists = mc.cmp(csm, csm, msr="MKL", prior=ap.prior)
pdists.tofile(f"pdists.mkl.prior{ap.prior}.f{suf}")
pdists = mc.cmp(csm, csm, msr="UWLLR", prior=ap.prior)
pdists.tofile(f"pdists.uwllr.prior{ap.prior}.f{suf}")
pdists = mc.cmp(csm, csm, msr="HELLINGER", prior=ap.prior)
pdists.tofile(f"pdists.hell.prior{ap.prior}.f{suf}")
#pdists = mc.cmp(csm, csm, msr="MKL", prior=0.1)
#pdists.tofile("pdists.mkl.prior0.1.f" + "32" if pdists.dtype.char == 'f' else "64")
#pdists = mc.cmp(csm, csm, msr="LLR", prior=0.1)
#pdists.tofile("pdists.llr.prior0.1.f" + "32" if pdists.dtype.char == 'f' else "64")
