#!/usr/bin/python
# -*- coding:utf-8 -*-

import os
import sys
import subprocess
import re

rel_dir = os.path.relpath(sys.argv[1])
bin_dir = os.path.abspath(rel_dir)
verbose = False
if len(sys.argv) > 2 and sys.argv[2] in ['-v', '--verbose']:
    verbose = True

for fname in os.listdir(bin_dir):
    exe = bin_dir + '/' + fname
    if os.path.isfile(exe) and os.access(exe, os.X_OK):
        title = fname.replace('_', ' ')            \
                     .replace('arma', 'Arma', 1)   \
                     .replace('eigen', 'Eigen', 1) \
                     .replace('vs', 'vs.', 1)
        title = re.sub(r'^([a-z]+)', lambda match: '['+match.group(1).upper()+']', title)
        cmd = [exe, '-r', 'html', '-o', exe + '.html', '-t', title]
        if verbose: cmd.append('-v')

        print('\n\033[0;32m' + '#'*len(title))
        print(title)
        print('#'*len(title) + '\033[0m\n')

        ret = subprocess.check_call(cmd)

        if verbose: print()

        print('\033[1;32m[report]  ====>  ' + rel_dir + '/' + fname + '.html\033[0m')

print('\n\033[1;32mall done\033[0m')
