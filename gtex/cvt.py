import numpy as np
import scipy.sparse as sp

# Script for converting back


ctsv = np.memmap("GTex.cts.u32", np.uint32)
is32 = True
if np.max(ctsv) < 65536:
    ctsv = ctsv.astype(np.uint16)
    is32 = False
mat = sp.csr_matrix(ctsv, np.memmap("GTex.ids.u16", np.uint16), np.memmap("GTex.indptr.u64", np.uint64))).T.tocsr()
mat.data = mat.data.astype(np.uint32 if is32 else np.uint16)
mat.data.tofile("GTex.T.cts.u" + "32" if is32 else "16")
mat.indices.astype(np.uint32).tofile("GTex.T.ids.u32")
mat.indptr.astype(np.uint64).tofile("GTex.T.indptr.u64")
np.array(mat.shape, np.uint64).tofile("GTex.T.shape.u64")
