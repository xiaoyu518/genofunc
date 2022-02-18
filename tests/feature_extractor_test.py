"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.feature_extractor import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestFeatureExtractor(unittest.TestCase):
    def test_run_feature_extractor(self):
        annotated_file = "%s/sequences/annotated.json" %data_dir
        gene = "gag pol"
        strip_gap = True
        filter_span = 0.9
        out_prefix = "%s/output/extracted_" %data_dir
        log_file = "%s/output/tmp.feature_extractor.log" %data_dir
        feature_extractor(annotated_file, gene, strip_gap, filter_span, out_prefix, log_file)
        os.unlink(log_file)
