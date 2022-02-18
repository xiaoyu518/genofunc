"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.filter_fasta import *

def run(options):

    filter_fasta(options.in_dir,
                 options.in_metadata,
                 options.gene_list,
                 options.min_length,
                 options.symmetric,
                 options.out_dir,
                 options.log_file
                 )