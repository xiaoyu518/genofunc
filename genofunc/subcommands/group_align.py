"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.group_align import *

def run(options):

    group_align(
        options.in_dir,
        options.group_size,
        options.reference_directory,
        options.out_dir,
        options.log_file
    )