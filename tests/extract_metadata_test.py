"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest

from genofunc.extract_metadata import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data')

class TestExtractMetadata(unittest.TestCase):
    def test_run_extract_metadata(self):
        in_fasta = "%s/sequences/seqB.fasta" %data_dir
        in_metadata = "%s/metadata/metadataB.csv" %data_dir
        column = ["country"]
        index_field = "strain"
        out_fasta = "%s/output/tmp.extract.fasta" %data_dir
        out_metadata = "%s/output/tmp.extracted_metadata.csv" %data_dir
        log_file = "%s/output/extract_metadata.log" %data_dir
        extract_metadata(in_fasta, in_metadata, column, index_field, out_fasta, out_metadata, log_file)
        os.unlink(out_fasta)
        os.unlink(out_metadata)
        os.unlink(log_file)