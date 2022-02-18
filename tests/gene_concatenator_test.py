"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.gene_concatenator import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestGeneConcatenator(unittest.TestCase):
    def test_run_gene_concatenator(self):
        in_fasta = ["%s/sequences/gene_fasta/pol.fasta", "%s/sequences/gene_fasta/gag.fasta" %data_dir]
        out_fasta = "%s/output/tmp.gag_pol.fasta" %data_dir
        log_file = "%s/output/mp.gene_concatenator.log" %data_dir
        gene_concatenator(in_fasta, out_fasta, log_file)
        os.unlink(out_fasta)
        os.unlink(log_file)