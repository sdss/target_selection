# encoding: utf-8
#
# @Author: Tom Dwelly
# @Date: Oct 2019
# @Filename: print_func.py
# @License: BSD 3-Clause
# @Copyright: Tom Dwelly

from __future__ import print_function

import sys
import os

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def print_error(*args):
    for ss in args:
        eprint('{} Error:   {} '.format(os.path.basename(sys.argv[0]),ss))

def print_warning(*args):
    for ss in args:
        eprint('{} Warning: {} '.format(os.path.basename(sys.argv[0]),ss))

def print_warning_stdout(*args):
    for ss in args:
        print('{} Warning: {} '.format(os.path.basename(sys.argv[0]),ss))

def print_comment(*args):
    for ss in args:
        print('{} Comment: {} '.format(os.path.basename(sys.argv[0]),ss))
