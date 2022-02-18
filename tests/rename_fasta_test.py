"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.rename_fasta import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestRenameFasta(unittest.TestCase):
    def test_run_rename_fasta(self):
        in_fasta = "%s/sequences/seq_pipe.fasta" %data_dir
        pipe = "|"
        out_fasta = "%s/output/depipe.fasta" %data_dir
        log_file = "%s/output/tmp.rename_fasta.log" %data_dir
        rename_fasta(in_fasta, pipe, out_fasta, log_file)
        os.unlink(out_fasta)
        os.unlink(log_file)

