import subprocess
from subprocess import check_output

def exact_wjaccard(p1, p2, k=17):
    return float(check_output("dashing2 sketch -k{k} --countdict --exact-kmer-dist --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])

def exact_jaccard(p1, p2, k=17):
    return float(check_output("dashing2 sketch -k{k} --set --exact-kmer-dist --cache --cmpout /dev/stdout {p1} {p2}").decode().strip('\n').split('\n')[-2].split("\t")[1])
