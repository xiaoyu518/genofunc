"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.genome_annotator import *

def run(options):

    genome_annotator(
        options.raw_fasta,
        options.reference_sequence,
        options.threads,
        options.out_annotation,
        options.log_file
    )
