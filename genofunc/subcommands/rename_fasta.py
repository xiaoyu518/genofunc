"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.rename_fasta import *

def run(options):

    rename_fasta(
        options.in_fasta,
        options.piping_character,
        options.out_fasta,
        options.log_file
    )