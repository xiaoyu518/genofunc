"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.reference_matcher import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestReferenceMatcher(unittest.TestCase):
    def test_run_reference_matcher(self):
        in_fasta = "%s/sequences/seqA.fasta" %data_dir
        reference_sequence = "%s/reference/reference.fasta" %data_dir
        out_fasta = "%s/output/referenced.fasta" %data_dir
        log_file = "%s/output/tmp.reference_matcher.log" %data_dir
        reference_matcher(in_fasta, reference_sequence, out_fasta, log_file)
        os.unlink(out_fasta)
        os.unlink(log_file)
