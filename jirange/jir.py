import subprocess
from subprocess import check_output

def exact_wjaccard(p1, p2, k=17):
    return float(check_output("dashing2 sketch -k{k} --countdict --exact-kmer-dist --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])

def exact_jaccard(p1, p2, k=17):
    return float(check_output("dashing2 sketch -k{k} --set --exact-kmer-dist --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])


def setsketch_jaccard(p1, p2, size, k=17, nb=8, fss=False):
    '''k is the k to use, and nb is the number of bytes to use in the setsketch-based set similarity estimation
       nb is 8 by default, which is Dashing 2's typical size
    '''
    return float(check_output("dashing2 sketch -k{k} {'--full-setsketch ' if fss else ''} -S{size} --fastcmp {nb} --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])
