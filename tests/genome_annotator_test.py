"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.genome_annotator import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestGeneAnnotator(unittest.TestCase):
    def test_run_genome_annotator(self):
        raw_fasta = "%s/sequences/referenced.fasta" %data_dir
        reference_sequence = "%s/reference/reference.json" %data_dir
        out_annotation = "%s/output/annotated.json" %data_dir
        log_file = "%s/output/tmp.genome_annotator.log" %data_dir
        genome_annotator(raw_fasta, reference_sequence, out_annotation, log_file)
        os.unlink(out_annotation)
        os.unlink(log_file)
