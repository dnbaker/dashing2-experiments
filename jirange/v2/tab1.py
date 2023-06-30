#!/usr/bin/env python3

import re

# from https://stackoverflow.com/questions/14693701/how-can-i-remove-the-ansi-escape-sequences-from-a-string-in-python

ansi_escape_8bit = re.compile(br'''
    (?: # either 7-bit C1, two bytes, ESC Fe (omitting CSI)
        \x1B
        [@-Z\\-_]
    |   # or a single 8-bit byte Fe (omitting CSI)
        [\x80-\x9A\x9C-\x9F]
    |   # or CSI + control codes
        (?: # 7-bit CSI, ESC [
            \x1B\[
        |   # 8-bit CSI, 9B
            \x9B
        )
        [0-?]*  # Parameter bytes
        [ -/]*  # Intermediate bytes
        [@-~]   # Final byte
    )
''', re.VERBOSE)

expected_columns = ['bin', 'k', 'sketchsize', 'n', 'Dash1', 'FSS1', 'FSS8', 'Mash', 'SS1', 'SS8']
use_sketches = ['8192', '32768', '131072']
use_tools = ['Mash', 'Dash1', 'FSS1', 'SS1']
use_ks = ['21', '31']

tools_rename = {'Dash1': 'Dashing1',
                'SS1': 'D2',
                'FSS1': 'D2-full',
                'SS8': 'D2-8byte'}

tools_renamed = []
for t in use_tools:
    if t in tools_rename:
        tools_renamed.append(tools_rename[t])
    else:
        tools_renamed.append(t)

sketches_rename = {'8192': '8',
                   '32768': '32',
                   '65536': '64',
                   '131072': '128',
                   '262144': '256',
                   '524288': '512'}

def find_in_string_list(ls, query):
    for i, item in enumerate(ls):
        if item == query:
            return i
    return None

print('\\definecolor{darkred}{rgb}{0.5,0,0}')
print('\\begin{table}')
print('\\begin{tabular}{rrrrrr}')
print(' & '.join(['k', 'kbits'] + tools_renamed) + ' \\\\ \\hline')
last_k = ''
with open('20sedf_ji_sse.tbl', 'rb') as fh:
    for ln in fh:
        ln = ln.rstrip()
        ln = ansi_escape_8bit.sub(b'', ln)
        if ln.startswith(b'#'):
            continue
        ln = ln.decode()
        toks = ln.strip().split()
        if toks[0] == 'bin':
            # first line
            if toks != expected_columns:
                raise ValueError('Expected ' + str(expected_columns) + ' got ' + str(toks))
            continue
        elif toks[0] == '<dbl>':
            # second line
            continue
        toks = toks[1:]
        k = toks[1]
        if k in use_ks and toks[2] in use_sketches:
            print_toks = []
            if k != last_k:
                print_toks.append(k)
                print_toks.append(toks[2])
            else:
                print_toks.append('')
                print_toks.append(toks[2])
            last_k = k
            results = []
            results_float = []
            for tool in use_tools:
                results.append(toks[find_in_string_list(expected_columns, tool)])
                results_float.append(float(results[-1]))
            min1 = list(sorted(results_float))[0]
            min2 = list(sorted(results_float))[1]
            for i in range(len(results)):
                if results_float[i] == min1:
                    results[i] = '{\color{red}' + results[i] + '}'
                elif results_float[i] == min2:
                    results[i] = '{\color{darkred}' + results[i] + '}'
            print_toks += results
            if print_toks[1] in sketches_rename:
                print_toks[1] = sketches_rename[print_toks[1]]
            print(' & '.join(print_toks) + ' \\\\')
print('\\end{tabular}')
print('\\caption{Sum of squared error JI}')
print('\\end{table}')
