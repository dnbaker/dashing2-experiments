import lzma
from glob import glob
import gzip
import numpy as np
from scipy.spatial.distance import squareform as SF

tissuesbin = np.fromfile("labels.tissue.u8.bin", 'B')
tissues = [x.split()[0] for x in open("labels.txt")]
tids = sorted(set(tissues))
NT = len(tids)
assert NT == 32
#NT = 1

def load_distmat(x):
    if x.endswith(".xz"):
        dat = lzma.open(x, "rb").read()
    elif x.endswith(".gz"):
        dat = gzip.open(x, "rb").read()
    else:
        dat = open(x, "rb").read()
    return SF(np.frombuffer(dat, np.float32))

def load_simmat(x):
    if x.endswith(".xz"):
        dat = lzma.open(x, "rb").read()
    elif x.endswith(".gz"):
        dat = gzip.open(x, "rb").read()
    else:
        dat = open(x, "rb").read()
    return SF(np.frombuffer(dat, np.uint16))

simmats = sorted(glob("*distmat.u[0-9][0-9]*xz"))
#simmats = ["teeny.distmat.u16"]
# assert len(simmats) == len(distmats), f"{len(simmats)}, {len(distmats)}"


# dms = list(map(load_distmat, distmats))
sms = list(map(load_simmat, simmats))
simnrgs = [int(x.split(".xz")[0].split(".")[-4]) for x in simmats]
tissue_masks = [np.where(tissuesbin == x) for x in np.arange(NT)]


def mean_noninf_distance(dmidx, i, j):
    lmask = tissue_masks[i]
    rmask = tissue_masks[j]
    sm = sms[dmidx]
    submat = sm[lmask][:,rmask[0]]
    filtered_submat = submat[submat > 0]
    dist_submat = -np.log(filtered_submat / simnrgs[dmidx])
    return np.mean(dist_submat)


def mean_sim(dmidx, i, j):
    lmask = tissue_masks[i]
    rmask = tissue_masks[j]
    sm = sms[dmidx]
    submat = sm[lmask][:,rmask[0]]
    return np.mean(submat)


def mean_nonzero_sim(dmidx, i, j):
    lmask = tissue_masks[i]
    rmask = tissue_masks[j]
    sm = sms[dmidx]
    submat = sm[lmask][:,rmask[0]]
    return np.mean(submat[submat > 0])

def make_xmat(dmidx, fn):
    return np.array([[fn(dmidx, i, j) for j in range(NT)] for i in range(NT)])


ninfdistmats = [make_xmat(x, mean_noninf_distance) for x in range(len(sms))]
nmatchmats = [make_xmat(x, mean_sim) for x in range(len(sms))]
nzeromatchmats = [make_xmat(x, mean_nonzero_sim) for x in range(len(sms))]

div = 1. / NT

for (i, path), ninfmat, nmatchmat, nzmatchmat in zip(enumerate(simmats), ninfdistmats, nmatchmats, nzeromatchmats):
    ncorrect_infdist = np.sum(np.argmin(ninfmat, axis=0) == np.arange(NT))
    ncorrect_match = np.sum(np.argmin(nmatchmat, axis=0) == np.arange(NT))
    ncorrect_zeromatch = np.sum(np.argmin(nzmatchmat, axis=0) == np.arange(NT))
    print(f"{path}\t{ncorrect_infdist}/{ncorrect_infdist * div}\t{ncorrect_match}/{ncorrect_match * div}\t{ncorrect_zeromatch}/{ncorrect_zeromatch * div}", flush=True)
