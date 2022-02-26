"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.merge import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestMerge(unittest.TestCase):
    def test_run_merge(self):
        in_fasta = ["%s/sequences/seqA.fasta", "%s/sequences/seqB.fasta" %data_dir]
        in_metadata = ["%s/metadata/metadataA.tsv", "%s/metadata/metadataB.tsv" %data_dir]
        index_field = "strain"
        out_dir = "%s/output/" %data_dir
        log_file = "%s/output/tmp.merge.log" %data_dir
        merge(in_fasta, in_metadata, index_field, out_dir, log_file)
        os.unlink(log_file)
