"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.reference_matcher import *

def run(options):

    reference_matcher(
        options.in_fasta,
        options.reference_sequence,
        options.out_fasta,
        options.log_file
    )