"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.concatenate_fasta import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestConcatenateFasta(unittest.TestCase):
    def test_run_concatenate_fasta(self):
        in_prefix = "%s/sequences/gene_fasta/" %data_dir
        gene = "pol gag"
        out_dir = "%s/output/" %data_dir
        log_file = "%s/output/tmp.concatenate_fasta.log" %data_dir
        concatenate_fasta(in_prefix, gene, out_dir, log_file)
        os.unlink(log_file)