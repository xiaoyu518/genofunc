"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.merge import *

def run(options):

    merge(
        options.in_fasta,
        options.in_metadata,
        options.index_field,
        options.out_dir,
        options.log_file
    )