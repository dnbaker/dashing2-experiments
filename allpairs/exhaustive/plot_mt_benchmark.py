import numpy as np
import sys
from itertools import zip_longest
import matplotlib.pyplot as plt

lines = list(filter(lambda x: x and x[0] != '#', map(str.strip, open("simple-filtered-stacked-subset.tsv"))))

toks = [x.split() for x in lines]; labels = [x[0] for x in toks];sizes = np.fromiter(map(int, (x[1] for x in toks)), np.int16);


labels = [x.replace("PMH", "D2W") for x in labels]
labels = np.array(labels)
nregs = np.fromiter(map(float, (x[2] for x in toks)), np.int16)
sketchbytes = np.fromiter(map(float, (x[3] for x in toks)), np.float16)

sketchtimes = np.fromiter(map(float, (x[5] for x in toks)), np.float32)

disttimes = np.fromiter(map(float, (x[6] for x in toks)), np.float32)

hasnobin = np.logical_not(np.fromiter(("bin" in x for x in labels), np.bool_))


CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
types = {"bindash", "mash", "d2fss", "d2op", "pmh", "d1"}
colors = {key: CB_color_cycle[x] for x, key in enumerate(types)}

def name2col(x):
    if "Bindash" in x:
        return "bindash"
    if "Dashing1" in x: return "d1"
    if "D2OP" in x: return "d2op"
    if "D2FSS" in x: return "d2fss"
    if "pmh" in x.lower() or "d2w" in x.lower(): return "pmh"
    return "mash"
    

ssz = sorted(set(nregs))
nbytes = sorted(set(sketchbytes))

rot = 45

for sz in ssz:
    szrows = np.where(np.logical_and(hasnobin, nregs == sz))[0]
    print(f"{len(szrows)} rows matching", file=sys.stderr)
    szlabels, szbytes, szdist, szst, szbytes = (x[szrows] for x in (labels, sketchbytes, disttimes, sketchtimes, sketchbytes))
    plt.clf()
    title = f"Parallel sketching benchmark - complete RefSeq, S = {sz}"
    sketchfn = f"cmpt-refseq-sketch.{sz}.all.svg"
    sketchorder = np.argsort(szst)
    nitems = len(szlabels)
    xticks = np.arange(nitems)
    fig, ax = plt.subplots()
    ax.set_xticks(xticks)
    szcolors = [colors[name2col(x)] for x in szlabels]
    assert len(szlabels) == len(xticks)
    ax.set_xticklabels(szlabels, rotation=rot, fontsize=4)
    bars = plt.bar(xticks, szst)
    for color, subbar in zip_longest(szcolors, bars):
        assert color is not None
        assert subbar is not None
        subbar.set_color(color)
    ax.set_xlabel("Sketching method")
    ax.set_ylabel("Sketch time, (s).")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(sketchfn)
    plt.savefig(sketchfn.replace("svg", "png"), dpi=300)
    szrows = np.where(szbytes == 8)

    plt.clf()
    fig, ax = plt.subplots()
    ax.set_xticks(xticks)
    szcolors = [colors[name2col(x)] for x in szlabels]
    ax.set_xticklabels(szlabels, rotation=rot, fontsize=4)
    bars = plt.bar(xticks, szdist)
    for color, subbar in zip_longest(szcolors, bars):
        assert color is not None
        assert subbar is not None
        subbar.set_color(color)
    ax.set_xlabel("Sketching method")
    ax.set_ylabel("Distance time, (s).")
    title = f"Parallel distance benchmark - complete RefSeq, S = {sz}"
    plt.title(title)
    sketchfn = sketchfn.replace("sketch", "cmp")
    plt.savefig(sketchfn)
    plt.savefig(sketchfn.replace("svg", "png"), dpi=300)

    plt.clf()
    fig, ax = plt.subplots()
    ax.set_xticks(xticks)
    szcolors = [colors[name2col(x)] for x in szlabels]
    ax.set_xticklabels(szlabels, rotation=rot, fontsize=4)
    bars = plt.bar(xticks, np.max(szdist) / szdist)
    for color, subbar in zip_longest(szcolors, bars):
        assert color is not None
        assert subbar is not None
        subbar.set_color(color)
    ax.set_xlabel("Sketching method")
    ax.set_ylabel("Speed-up time, relative to slower (usu. Mash)")
    title = f"Parallel distance benchmark - complete RefSeq, S = {sz}"
    plt.title(title)
    sketchfn = f"cmpt-refseq-cmp-relative.{sz}.all.svg"
    plt.savefig(sketchfn)
    plt.savefig(sketchfn.replace("svg", "png"), dpi=300)


    szlabels, szbytes, szdist, szst, szbytes = (x[szrows] for x in (szlabels, szbytes, szdist, szst, szbytes))
    plt.clf()
    title = f"Parallel sketching - complete RefSeq, S = {sz}"
    sketchfn = f"cmpt-refseq-sketch.{sz}.fullregsonly.svg"
    sketchorder = np.argsort(szst)
    nitems = len(szlabels)
    xticks = np.arange(nitems)
    fig, ax = plt.subplots()
    ax.set_xticks(xticks)
    szcolors = [colors[name2col(x)] for x in szlabels]
    assert len(szlabels) == len(xticks)
    ax.set_xticklabels(szlabels, rotation=rot, fontsize=4)
    bars = plt.bar(xticks, szst)
    for color, subbar in zip_longest(szcolors, bars):
        assert color is not None
        assert subbar is not None
        subbar.set_color(color)
    ax.set_xlabel("Sketching method")
    ax.set_ylabel("Sketch time, (s).")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(sketchfn)
    plt.savefig(sketchfn.replace("svg", "png"), dpi=300)


    plt.clf()
    fig, ax = plt.subplots()
    ax.set_xticks(xticks)
    szcolors = [colors[name2col(x)] for x in szlabels]
    ax.set_xticklabels(szlabels, rotation=rot, fontsize=4)
    bars = plt.bar(xticks, szdist)
    for color, subbar in zip_longest(szcolors, bars):
        assert color is not None
        assert subbar is not None
        subbar.set_color(color)
    ax.set_xlabel("Sketching method")
    ax.set_ylabel("Distance time, (s).")
    title = f"Parallel distance benchmark - complete RefSeq, S = {sz}"
    plt.title(title)
    sketchfn = f"cmpt-refseq-cmp.{sz}.fullregsonly.svg"
    plt.savefig(sketchfn)
    plt.savefig(sketchfn.replace("svg", "png"), dpi=300)


