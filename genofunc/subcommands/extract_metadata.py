"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
"""

from genofunc.extract_metadata import *

def run(options):

    extract_metadata(options.in_fasta,
                     options.in_metadata,
                     options.column_name,
                     options.column_rename,
                     options.id_column,
                     options.out_fasta,
                     options.out_metadata,
                     options.log_file
                     )