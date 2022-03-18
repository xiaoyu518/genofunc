"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.feature_extractor import *

def run(options):

    feature_extractor(
        options.in_annotation,
        options.gene_region,
        options.strip_gap,
        options.filter_coverage,
        options.filter_span,
        options.out_prefix,
        options.log_file
    )
