"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from genofunc.strain_encoder import *


def run(options):

    strain_encoder(
        options.in_fasta,
        options.in_metadata,
        options.encoding_column,
        options.out_dir,
        options.log_file
    )
