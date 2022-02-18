"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.strain_encoder import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestStrainEncoder(unittest.TestCase):
    def test_run_strain_encoder(self):
        in_fasta = "%s/sequences/seqB.fasta" %data_dir
        in_metadata = "%s/metadata/seq_pipe.fasta" %data_dir
        encoding_column = 2
        out_dir = "%s/output/" %data_dir
        log_file = "%s/output/tmp.strain_encoder.log" %data_dir
        strain_encoder(in_fasta, in_metadata, encoding_column, out_dir, log_file)
        os.unlink(log_file)
