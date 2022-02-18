"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.name_splitter import *

def run(options):

    name_splitter(options.in_fasta,
                  options.piping_character,
                  options.header,
                  options.out_metadata,
                  options.log_file)
