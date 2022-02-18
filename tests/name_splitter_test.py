"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.name_splitter import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestNameSplitter(unittest.TestCase):
    def test_run_name_splitter(self):
        in_fasta = "%s/sequences/seq_pipe.fasta" %data_dir
        pipe = "|"
        out_metadata = "%s/output/tmp.output_metadata.csv" %data_dir
        log_file = "%s/output/tmp.name_splitter.log" %data_dir
        name_splitter(in_fasta, pipe, out_metadata, log_file)
        os.unlink(out_metadata)
        os.unlink(log_file)