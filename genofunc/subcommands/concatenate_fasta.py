"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.concatenate_fasta import *

def run(options):

    concatenate_fasta(options.in_prefix,
                      options.gene_region,
                      options.out_dir,
                      options.log_file)
