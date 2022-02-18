"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.filter_fasta import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestFilterFasta(unittest.TestCase):
    def test_run_filter_fasta(self):
        in_dir = "%s/sequences/gene_fasta" %data_dir
        in_metadata = "%s/metadata/metadataB.csv" %data_dir
        gene = "pol gag"
        min_length = 0.9
        symmetric = True
        out_dir = "%s/output/" %data_dir
        log_file = "%s/output/tmp.filter_fasta.log" %data_dir
        filter_fasta(in_dir, in_metadata, gene, min_length, symmetric, out_dir, log_file)
        os.unlink(log_file)
