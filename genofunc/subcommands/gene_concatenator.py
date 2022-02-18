"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.gene_concatenator import *

def run(options):

    gene_concatenator(options.in_fasta,
                      options.out_fasta,
                      options.log_file)
