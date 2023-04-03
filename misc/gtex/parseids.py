'''
This extracts sample IDs from the GTex experiment, using
recountcsr's output parsed.remap and the samples.tsv file
from Recount.
'''


def tok52sex(x):
    if not x:
        return 0
    return int(x)

def parseline(x):
    if x.startswith("rail_id"):
        return tuple()
    '''Returns railid, exid, tissue type, sex, and age'''
    toks = x.split('\t')
    return (int(toks[0]), toks[2], toks[3], int(toks[5]) if toks[5] else 0, toks[6])

def tups2id(tup):
    return ":".join(map(str, (tup[2], tup[1], tup[3], tup[0], tup[3], tup[4])))

remapdict, invmapdict = {}, {}
for line in open("parsed.remap"):
    rail, distid = map(int, line.strip().split(':'))
    remapdict[rail] = distid
    invmapdict[distid] = rail

lines = list(filter(lambda x: x, map(parseline, open("samples.tsv"))))
res = ["" for i in range(19214)]
for l in lines:
    rid, exid, tissue, sex, age = l
    res[remapdict[rid]] = tups2id(l)

assert all(res)

for r in res:
    print(r)
