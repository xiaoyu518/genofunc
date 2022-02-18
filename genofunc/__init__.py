"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2022 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk).
"""

from pkg_resources import get_distribution

try:
    __version__ = get_distribution("genofunc").version
except:
    __version__ = "local"

__all__ = ["concatenate_fasta", "extract_metadata", "feature_extractor", "filter_fasta", "gene_concatenator",
           "genome_annotator", "group_align","merge","name_splitter", "reference_matcher", "rename_fasta", 
           "strain_encoder"]

from genofunc import *
