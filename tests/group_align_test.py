"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.group_align import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestGroupAlign(unittest.TestCase):
    def test_run_group_align(self):
        in_dir = "%s/sequences/group_align/" %data_dir
        group_size = 10
        reference_dir = "%s/reference/" %data_dir
        out_dir = "%s/output/" %data_dir
        log_file = "%s/output/tmp.group_align.log" %data_dir
        group_align(in_dir, group_size, reference_dir, out_dir, log_file)
        os.unlink(log_file)
